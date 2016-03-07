% McDermott
% 10-14-2015
% poiscell_2d_init.m

% First, create a 3 x 3 array of cells with a nontrivial solenoidal velocity field.

Lx = 2*pi;
Ly = 2*pi;
Sx = rand(1)*Lx; % random phase shift in x
Sy = rand(1)*Ly; % random phase shift in y

nx = 3;
ny = 3;
B = 2;

% mesh spacing
dx = Lx/nx;
dy = Ly/ny;

% staggered face locations
x = [0:nx]*dx + Sx;
y = [0:ny]*dy + Sy;

% cell center locations
xp = x(1:nx) + 0.5*dx;
yp = y(1:ny) + 0.5*dy;

for i=1:nx+1
    for j=1:ny
        uhat(i,j) = 1 - B*cos(x(i))*sin(yp(j));
    end
end

for i=1:nx
    for j=1:ny+1
        vhat(i,j) = 1 + B*sin(xp(i))*cos(y(j));
    end
end

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

% solve Poisson equation for "condensed" system
n_cells = nx*ny;
Ac = A(1:n_cells-1,1:n_cells-1);
bc = b(1:n_cells-1);
pc = Ac\bc';
pvec = [pc;0];
    
% map solution vector to computational indices
for i=1:nx
    for j=1:ny
        np = (j-1)*nx + i;
        p(i,j) = pvec(np);
    end
end

% project velocities
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

for j=1:ny
    u(nx+1,j) = u(1,j);
end

for i=1:nx
    v(i,ny+1) = v(i,1);
end








