% McDermott
% 8-26-2021
% channel_cc_ibm.m
%
% Simple staggered primitive variable treatment channel flow
% using a simple FE time integrator.

close all
clear all

Lx = 2; % channel dimension in x
Ly = 1; % channel dimension in y
nx = 16; % number of pcells in x
ny = 8; % number of pcells in y
dx = Lx/nx; % uniform grid spacing in x
dy = Ly/ny; % uniform grid spacing in y
dxdx = dx^2;
dydy = dy^2;
U_INLET = 1; % inlet velocity
nu = 1e-5; % kinematic viscosity (1/Re)
Fo = 0.5; % Fourier number
CFL = 0.5; % Courant number
T = 10; % total simulation time

% For the CC_IBM test, let's blank out cells on the top and bottom
% of the channel.

ncc = 2; % number of completely solid cells on top and bottom of channel
cfa = 1; % cutface area factor

% STAGGERED GRID ARRANGEMENT:
%
% Let this represent the bottom left pcell control volume.
% Velocities are stored on their respective faces.  Normal stresses
% and normal advective fluxes are stored at the pcell center.
% Off-diagonal stresses and advective fluxes are stored at vertices
% marked by the Xs.  The bottom left most element is prescribed the
% indices (1,1).  Because of this, some of the differencing
% below looks confusing, but this seems somewhat unavoidable.
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
%       |   FZX(1,1)       H(1,1)                |   FZX(2,1)
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

% initialization
u = zeros(nx+1,ny);
v = zeros(nx,ny+1);
FZX = u;
FZY = v;
u_hat = zeros(nx+1,ny);
v_hat = zeros(nx,ny+1);
up = zeros(nx,ny);
vp = zeros(nx,ny);
Z = zeros(nx,ny);
ZPAD = zeros(nx+1,ny+1); % pad for pcolor output
kres = zeros(nx,ny);
tau11 = zeros(nx,ny);
tau22 = zeros(nx,ny);
tau12 = zeros(nx+1,ny+1);
omega3 = zeros(nx+1,ny+1);
H = zeros(nx,ny);
b_vec = zeros(nx*ny,1);
HNXP1 = zeros(1,ny);

% cell-centered grid
x = linspace(dx/2,Lx-dx/2,nx);
y = linspace(dy/2,Ly-dy/2,ny);
for j = 1:ny
    X(:,j) = x;
end
for i = 1:nx
    Y(i,:) = y;
end

celltype = ones(nx+1,ny+1); % pad this array for later use with pcolor output
facextype = ones(nx+1,ny);
faceytype = ones(nx,ny+1);
if ncc>0
    celltype(:,1:ncc) = 0;
    celltype(:,(ny-ncc+1):ny+1) = 0;

    % cutcells
    celltype(:,ncc+1) = cfa;
    celltype(:,ny-ncc) = cfa;
end
for i=1:nx
    for j=1:ny
        if celltype(i,j)==0
            facextype(i,j) = 0; facextype(i+1,j) = 0;
            faceytype(i,j) = 0; faceytype(i,j+1) = 0;
        end
        if celltype(i,j)>0 & celltype(i,j)<1
            facextype(i,j) = cfa; facextype(i+1,j) = cfa;  % cutface in x direction
            if j<floor(ny/2)
                faceytype(i,j) = 0; faceytype(i,j+1) = 1;  % bottom of cell is solid, top is a regular gas face
            else
                faceytype(i,j) = 1; faceytype(i,j+1) = 0;  % bottom of cell is gas, top solid
            end
        end
    end
end

% set inlet bc
U0 = zeros(1,ny);
Z0 = ones(1,ny);
for j=1:ny
    U0(j)=U_INLET*facextype(1,j);
    u(1,j)=U0(j);
    if facextype(1,j)==0
        Z0(j)=0;
    end
end

% node-centered grid
xf = linspace(0,Lx,nx+1);
yf = linspace(0,Ly,ny+1);
for j = 1:ny+1
    XF(:,j) = xf;
end
for i = 1:nx+1
    YF(i,:) = yf;
end

% Build A
A = build_sparse_matrix_mixed([nx ny],[dx dy],[0 2 0 0]);

% Main time step loop
t = 0;
%T = 1e-6;
while t<T

    dt_dif = Fo/(nu*(1/dxdx+1/(dydy*cfa^2))) % time step based on Fourier number
    velmax = max([max(abs(U0)),max(max(abs(u))),max(max(abs(v)))])
    dt_adv = CFL*min(dx,dy*cfa)/velmax % time step based on CFL
    dt = min([dt_dif,dt_adv]) % time step used in Forward Euler integrator
    t = t+dt

    % scalar fluxes (Godunov)
    FZX(1,:)=Z0.*u(1,:);
    for j=1:ny
        for i=2:nx+1
            if u(i,j)>=0
                FZX(i,j)=Z(i-1,j)*u(i,j);
            else
                if i<=nx
                    FZX(i,j)=Z(i,j)*u(i,j);
                else
                    FZX(i,j)=0; % assume no tracer from outlet side
                end
            end
        end
    end

    FZY(:,1)=0;
    FZY(:,ny+1)=0;
    for j=2:ny
        for i=1:nx
            if v(i,j)>=0
                FZY(i,j)=Z(i,j-1)*v(i,j);
            else
                FZY(i,j)=Z(i,j)*v(i,j);
            end
        end
    end

    % update scalars
    for j=1:ny
        for i=1:nx
            Z(i,j) = Z(i,j) - dt*( (FZX(i+1,j)-FZX(i,j))/dx + (FZY(i,j+1)-FZY(i,j))/dy );
        end
    end

    % compute normal stresses
    for j = 1:ny
        for i = 1:nx
            tau11(i,j) = -2*nu*( u(i+1,j)-u(i,j) )/dx;
            tau22(i,j) = -2*nu*( v(i,j+1)-v(i,j) )/dy;

            kres(i,j) = 0.5*(up(i,j)^2 + vp(i,j)^2);
        end
    end

    % compute interior vorticity and off-diagonal stresses
    for j = 2:ny
        for i = 2:nx
            omega3(i,j) = (v(i,j)-v(i-1,j))/dx - (u(i,j)-u(i,j-1))/dy;
            tau12(i,j) = -nu*( (u(i,j)-u(i,j-1))/dy + (v(i,j)-v(i-1,j))/dx );
        end
    end

    % left wall stress and vorticity
    i = 1;
    for j = 2:ny
        dvdx_wall = 0;
        tau12(i,j) = -nu*dvdx_wall;
        omega3(i,j) = dvdx_wall;
    end
    % right wall stress and vorticity
    i = nx+1;
    for j = 2:ny
        dvdx_wall = 0;
        tau12(i,j) = -nu*dvdx_wall;
        omega3(i,j) = dvdx_wall;
    end

    % bottom wall stress and vorticity
    j = 1;
    for i = 1:nx+1
        dudy_wall = 2*u(i,j)/dy;
        tau12(i,j) = -nu*dudy_wall;
        omega3(i,j) = -dudy_wall;
    end
    % top wall stress and vorticity
    j = ny+1;
    for i = 1:nx
        dudy_wall = -2*u(i,j-1)/dy;
        tau12(i,j) = -nu*dudy_wall;
        omega3(i,j) = -dudy_wall;
    end

    % u momentum -- predictor
    for j = 1:ny
        for i = 2:nx

            Fu_visc = -( (tau11(i,j)-tau11(i-1,j))/dx + (tau12(i,j+1)-tau12(i,j))/dy );

            vn = 0.5*( v(i-1,j+1) + v(i,j+1) );
            vs = 0.5*( v(i-1,j)   + v(i,j)   );

            Fx = -0.5*( vs*omega3(i,j) + vn*omega3(i,j+1) );

            % immersed boundary forcing
            if facextype(i,j)==0
                Fu_visc = 0;
                Fx = (0 - u(i,j))/dt - ( H(i,j)-H(i-1,j) )/dx;
            end

            u_hat(i,j) = u(i,j) + dt*(Fx + Fu_visc);
        end

        % inlet boundary forcing
        Fx = (U0(j) - u(1,j))/dt;
        u_hat(1,j) = u(1,j) + dt*Fx;

        % outlet boundary forcing
        if facextype(i,j)==0
            Fu_visc = 0;
            Fx = (0 - u(nx+1,j))/dt - ( HNXP1(j)-H(nx,j) )/dx;
        else
            Fu_visc = 0;
            vn = v(nx,j+1);
            vs = v(nx,j);
            Fx = -0.5*( vs*omega3(nx,j) + vn*omega3(nx,j+1) );
        end
        u_hat(nx+1,j) = u(nx+1,j) + dt*(Fx + Fu_visc);
    end

    % v momentum -- predictor
    for j = 2:ny
        for i = 1:nx

            Fv_visc = -( (tau12(i+1,j)-tau12(i,j))/dx + (tau22(i,j)-tau22(i,j-1))/dy );

            ue = 0.5*( u(i+1,j-1) + u(i+1,j) );
            uw = 0.5*( u(i,j-1)   + u(i,j)   );

            Fy = 0.5*( uw*omega3(i,j) + ue*omega3(i+1,j) );

            % immersed boundary forcing
            if faceytype(i,j)==0
                Fv_visc = 0;
                Fy = (0 - v(i,j))/dt - ( H(i,j)-H(i,j-1) )/dy;
            end

            v_hat(i,j) = v(i,j) + dt*(Fy + Fv_visc);
        end
    end

    % build source (Poisson right-hand-side)
    for j = 1:ny
        for i = 1:nx
            b(i,j) = ( ( u_hat(i+1,j)-u_hat(i,j) )/dx + ( v_hat(i,j+1)-v_hat(i,j) )/dy ) / dt;
        end

        % apply Dirichlet bcs to outflow
        if u(nx+1,j)>0 & facextype(nx+1,j)>0
            b(nx,j) = b(nx,j) - 2*kres(nx,j)/dxdx;
        end
    end

    % map b to source vector
    for j = 1:ny
        for i = 1:nx
            p = (j-1)*nx+i;
            b_vec(p) = b(i,j);
        end
    end
    % % subtract mean for discrete compatibility condition
    % b_vec = b_vec - mean(b_vec);

    % solve for H vector
    H_vec = A\b_vec;
    % % subtract mean for arbitrary solution
    % H_vec = H_vec - mean(H_vec);
    H_min = min(H_vec);
    H_max = max(H_vec);

    % map H_vec to grid
    for j = 1:ny
        for i = 1:nx
            p = (j-1)*nx+i;
            H(i,j) = H_vec(p);
        end
    end

    % project velocities
    for j = 1:ny
        for i = 2:nx
            u(i,j) = u_hat(i,j) - dt*( H(i,j)-H(i-1,j) )/dx;
        end
        % apply inflow bc
        u(1,j) = u_hat(1,j);
        % apply outflow bc
        if u(nx+1,j)>0 & facextype(nx+1,j)>0
            HNXP1(j) = 2*kres(nx,j) - H(nx,j);
        else
            HNXP1(j) = -H(nx,j);
        end
        u(nx+1,j) = u_hat(nx+1,j) - dt*( HNXP1(j)-H(nx,j) )/dx;
    end
    for j = 2:ny
        for i = 1:nx
            v(i,j) = v_hat(i,j) - dt*( H(i,j)-H(i,j-1) )/dy;
        end
    end

    % check divergence
    for j = 1:ny
        for i = 1:nx
            b(i,j) = ( u(i+1,j)-u(i,j) )/dx + ( v(i,j+1)-v(i,j) )/dy;
        end
    end
    % u
    % b
    % return
    display(['max stage 1 velocity divergence = ',num2str( max(max(abs(b))) )])
    
    % interpolate velocities to cell centers
    for j = 1:ny
        for i = 1:nx
            up(i,j) = 0.5*( u(i+1,j)+u(i,j) );
            vp(i,j) = 0.5*( v(i,j+1)+v(i,j) );
        end
    end

    subplot(3,1,1)
    % pcolor(XF,YF,celltype)
    ZPAD(1:nx,1:ny) = Z;
    pcolor(XF,YF,ZPAD)
    colorbar
    title('Cutcell Mesh')
    axis([0 Lx 0 Ly])
    set(gca,'PlotBoxAspectRatio',[Lx Ly 1])
    
    subplot(3,1,2)
    velmag = sqrt( up.*up + vp.*vp );
    V = linspace(0,1,20);
    contourf(X,Y,velmag,V); hold on
    colorbar
    quiver(X,Y,up,vp); hold off
    title('velocity vectors')
    axis([0 Lx 0 Ly])
    set(gca,'PlotBoxAspectRatio',[Lx Ly 1])
    
    subplot(3,1,3)
    V = linspace(H_min,H_max,20);
    contourf(X,Y,H,V)
    title('pressure contours')
    colorbar
    set(gca,'PlotBoxAspectRatio',[Lx Ly 1])
    
    pause(0.001)

end

