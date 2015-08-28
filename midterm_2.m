% McDermott
% 3-3-13
% midterm_2.m

close all
clear all
warning('off','all');

% select formulation
% 1 = advective form
% 2 = conservative form
% 3 = Stokes form

formulation = 3;

% physcial parameters

B = 2;       % amplitude, m/s    
nu = 0.;     % kinematic viscosity, m2/s

Lx = 2*pi;   % x domain length, m
Ly = 2*pi;   % y domain length, m
nx = 32;     % number of cells in x
ny = 32;     % number of cells in y
T  = 2*pi;   % total time, s
cfl = 0.25;  % CFL

% mesh spacing
dx = Lx/nx;
dy = Ly/ny;

% staggered face locations
x = [0:nx]*dx;
y = [0:ny]*dy;

% cell center locations
xp = x(1:nx)+0.5*dx;
yp = y(1:ny)+0.5*dy;

% STAGGERED GRID ARRANGEMENT:
%
% Let this represent the bottom left pcell control volume.
% Velocities are stored on their respective faces.  Normal stresses
% and normal advective fluxes are stored at the pcell center.
% Off-diagonal stresses and advective fluxes are stored at vertices
% marked by the Xs.  The bottom left most element is prescribed the
% indices (1,1).
%
%       omega3(1,2)                              omega3(2,2)
%       tau12(1,2)          ^ v(1,2)             tau12(2,2)
%       X-------------------|--------------------X
%       |                                        |
%       |                                        |
%       |                                        |
%       |                                        |
%       |                                        |
%       |                                        |
%       |                                        |
%       |                                        |
%      ---> u(1,1)          O                   ---> u(2,1)
%       |                  p(1,1) or H(1,1)      |
%       |                  tau11(1,1)            |
%       |                  tau22(1,1)            |
%       |                                        |
%       |                                        |
%       |                                        |
%       |                                        |
%       |                   ^ v(1,1)             |
%       X-------------------|--------------------X
%       omega3(1,1)                              omega3(2,1)
%       tau12(1,1)                               tau12(2,1)
%
% Notes on the pressure term:
%
% In a convenctional formulation (advective, conservative), the pressure is
% denoted p and has the physical interpretation of being the isotropic part
% of the total stress (which we call the "mechanical pressure").  It is an
% assumption (called "Stokes' assumption") that the mechanical pressure is
% equivalent to the thermodynamic pressure (from, say, the ideal gas law).
%
% When the Navier-Stokes equations are written in Stokes form, the gradient
% of the kinetic energy per unit mass (k) appears and is combined with the
% mechanical pressure to create a stagnation energy per unit mass.  This is
% denoted H (by FDS) and is often cavalierly referred to as the
% "pseudo-pressure".

% initial condition

for i=1:nx
    for j=1:ny
        u(i,j) = 1 - B*cos(x(i))*sin(yp(j));
        v(i,j) = 1 + B*sin(xp(i))*cos(y(j));
    end
end

t = 0;
dt = cfl*dx/(B+1);

set(gcf,'DefaultAxesFontSize',16)
set(gcf,'DefaultTextFontSize',16)

% build A matrix (periodic)

A = sparse(nx*ny,nx*ny);
for i=1:nx
    for j=1:ny
        
        ip1=i+1;
        im1=i-1;
        jp1=j+1;
        jm1=j-1;
        
        if ip1>nx; ip1=ip1-nx; end
        if jp1>ny; jp1=jp1-ny; end
        if im1<1; im1=im1+nx; end
        if jm1<1; jm1=jm1+ny; end
        
        % lexicographical ordering
        np    = (j-1)*nx + i;
        east  = (j-1)*nx + ip1;
        west  = (j-1)*nx + im1;
        north = (jp1-1)*nx + i;
        south = (jm1-1)*nx + i;
        
        A(np,np   ) = -(2/dx^2 + 2/dy^2);
        A(np,east ) = 1/dx^2;
        A(np,west ) = 1/dx^2;
        A(np,north) = 1/dy^2;
        A(np,south) = 1/dy^2;
    end
end

% uncomment to view matrix structure

% full(A)
% spy(A)

% return
% tightfig(gcf);
% print(gcf,'-dpdf','spyA')

while t<T
    
    t = t + dt;
    
    for i=1:nx
        for j=1:ny
            
            ip1=i+1;
            im1=i-1;
            jp1=j+1;
            jm1=j-1;
        
            if ip1>nx; ip1=ip1-nx; end
            if jp1>ny; jp1=jp1-ny; end
            if im1<1;  im1=im1+nx; end
            if jm1<1;  jm1=jm1+ny; end
             
            % compute advective stress
            
            if formulation==1
                
                % advective form
                ududx(i,j) = 0.5*(u(i,j)+u(ip1,j))*(u(ip1,j)-u(i,j))/dx;
                vdvdy(i,j) = 0.5*(v(i,j)+v(i,jp1))*(v(i,jp1)-v(i,j))/dy;
                
                vdudy(i,j) = 0.5*(v(i,j)+v(im1,j))*(u(i,j)-u(i,jm1))/dy;
      
                udvdx(i,j) = 0.5*(u(i,j)+u(i,jm1))*(v(i,j)-v(im1,j))/dx;          
                
            elseif formulation==2
                
                % conservative form
                uu(i,j) = ( 0.5*(u(i,j)+u(ip1,j)) )^2; % cell center
                vv(i,j) = ( 0.5*(v(i,j)+v(i,jp1)) )^2; % cell center
                uv(i,j) = 0.5*(u(i,j)+u(i,jm1)) * 0.5*(v(i,j)+v(im1,j)); % vertex
                
            elseif formulation==3
                
                % Stokes form
                omega3(i,j) = (v(i,j)-v(im1,j))/dx - (u(i,j)-u(i,jm1))/dy;
                ubar(i,j) = 0.5*(u(i,jm1)+u(i,j));
                vbar(i,j) = 0.5*(v(im1,j)+v(i,j));
                
            end
            
            % compute viscous stress components
            dudy = (u(i,j)-u(i,jm1))/dy;
            dvdx = (v(i,j)-v(im1,j))/dx;
            dudx = (u(ip1,j)-u(i,j))/dx;
            dvdy = (v(i,jp1)-v(i,j))/dy;
            
            tau11(i,j) = -2*nu*dudx; % cell center
            tau22(i,j) = -2*nu*dvdy; % cell center
            tau12(i,j) = -nu*(dudy + dvdx); % vertex
            
        end
    end
    
    % compute force terms and update velocity predictor
    for i=1:nx
        for j=1:ny
            
            ip1=i+1;
            im1=i-1;
            jp1=j+1;
            jm1=j-1;
            
            if ip1>nx; ip1=ip1-nx; end
            if jp1>ny; jp1=jp1-ny; end
            if im1<1;  im1=im1+nx; end
            if jm1<1;  jm1=jm1+ny; end
            
            if formulation==1
                
                % advective form
                Fx = 0.5*( ududx(i,j)+ududx(im1,j) + vdudy(i,j)+vdudy(i,jp1) );
                Fy = 0.5*( udvdx(i,j)+udvdx(ip1,j) + vdvdy(i,j)+vdvdy(i,jm1) );
                
            elseif formulation==2
                
                % conservative form
                Fx = (uu(i,j)-uu(im1,j))/dx + (uv(i,jp1)-uv(i,j))/dy;
                Fy = (uv(ip1,j)-uv(i,j))/dx + (vv(i,j)-vv(i,jm1))/dy;
                
            elseif formulation==3
                
                % Stokes form
                Fx = -0.5*( vbar(i,j)*omega3(i,j) + vbar(i,jp1)*omega3(i,jp1) );
                Fy =  0.5*( ubar(i,j)*omega3(i,j) + ubar(ip1,j)*omega3(ip1,j) );
                
            end
            
            % add viscous stresses
            Fx = Fx + (tau11(i,j)-tau11(im1,j))/dx + (tau12(i,jp1)-tau12(i,j))/dy;
            Fy = Fy + (tau12(ip1,j)-tau12(i,j))/dx + (tau22(i,j)-tau22(i,jm1))/dy;
            
            uhat(i,j) = u(i,j) - dt*Fx;
            vhat(i,j) = v(i,j) - dt*Fy;
        end
    end
    
    % build right hand side of Poisson equation
    for i=1:nx
        for j=1:ny

            ip1=i+1;
            jp1=j+1;
            if ip1>nx; ip1=ip1-nx; end
            if jp1>ny; jp1=jp1-ny; end
            
            np = (j-1)*nx + i;
            b(np) = (uhat(ip1,j)-uhat(i,j))/dx + (vhat(i,jp1)-vhat(i,j))/dy;
        end
    end
    b = b-mean(b); % for discrete compatibility, b should have zero mean
    
    % solve Poisson equation
    pvec = A\b';
    
    % map solution vector to computational indices
    for i=1:nx
        for j=1:ny
            np = (j-1)*nx + i;
            p(i,j) = pvec(np);
        end
    end
    
    % project velocities
    % note: dt may be omitted here if it is left out of the b vector (right hand side) above
    for i=1:nx
        for j=1:ny
            
            im1=i-1;
            jm1=j-1;
            if im1<1; im1=im1+nx; end
            if jm1<1; jm1=jm1+ny; end
            
            u(i,j) = uhat(i,j) - ( p(i,j) - p(im1,j) )/dx;
            v(i,j) = vhat(i,j) - ( p(i,j) - p(i,jm1) )/dy;
            
        end
    end
    
    % check divergence
    
    for i=1:nx
        for j=1:ny

            ip1=i+1;
            jp1=j+1;
            if ip1>nx; ip1=ip1-nx; end
            if jp1>ny; jp1=jp1-ny; end
            
            div(i,j) = (u(ip1,j)-u(i,j))/dx + (v(i,jp1)-v(i,j))/dy;
        end
    end
    
    display( ['max divergence = ',num2str(max(max(abs(div))))] )
    
    surf(xp,yp,u)
    xlabel('x')
    ylabel('y')
    zlabel('u') 
    axis([0 Lx 0 Ly -1 3])
    
    pause(0.001)
    
end

% colorbar('Ylim',[-1,3])
% tightfig(gcf);
% print(gcf,'-dpdf','u_num')

% compute the error

for i=1:nx
    for j=1:ny
        
        u_exact = 1 - B*cos(x(i)-t)*sin(yp(j)-t)*exp(-2*nu*t);
        v_exact = 1 + B*sin(xp(i)-t)*cos(y(j)-t)*exp(-2*nu*t);     
        p_exact = -(B^2)/4*( cos(2*(xp(i)-t)) + cos(2*(yp(j)-t)) )*exp(-4*nu*t);
        
        if formulation==3
            up_exact = 1 - B*cos(xp(i)-t)*sin(yp(j)-t)*exp(-2*nu*t);
            vp_exact = 1 + B*sin(xp(i)-t)*cos(yp(j)-t)*exp(-2*nu*t);  
            p_exact = p_exact + 0.5*(up_exact^2 + vp_exact^2); % definition of H
        end
        
        uerr(i,j) = u(i,j) - u_exact;
        verr(i,j) = v(i,j) - v_exact;
        perr(i,j) = p(i,j) - p_exact;

    end
end

% subtract mean from pressure error

perr = perr - mean(mean(perr));

format long e
L2_uerr = norm(uerr)/(nx*ny)
L2_verr = norm(verr)/(nx*ny)
L2_perr = norm(perr)/(nx*ny)

figure
set(gcf,'DefaultAxesFontSize',16)
set(gcf,'DefaultTextFontSize',16)
surf(x(1:nx),yp,uerr)          
xlabel('x')
ylabel('y')
zlabel('uerr') 








