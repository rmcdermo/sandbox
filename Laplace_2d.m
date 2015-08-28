% McDermott
% 10-29-2014
% Laplace_2d.m

close all
clear all

Lx = 1;
Ly = 1;
nx = 8;
ny = 8;
dx = Lx/nx;
dy = Ly/ny;
x = dx/2:dx:Lx-dx/2;
y = dy/2:dy:Ly-dy/2;

% lexicographic mapping
for i=1:nx
    for j=1:ny
        p(i,j) = (j-1)*nx+i;
    end
end

A = build_sparse_matrix_2([nx ny],[dx dy],[0 2 0 0]);

% source term

b_vec = zeros(nx*ny,1);

% left boundary
s = 0.5*y;
i = 1;
for j = 1:ny
    b_vec(p(i,j)) = s(j)/dx;
end

% define block

nxb = (round(.375*nx)+1):(.625*nx);
nyb = nxb;

for i=nxb
    for j=nyb
        % modify i-1,i+1 coefficients
        A(p(i-1,j),p(i-1,j))=A(p(i-1,j),p(i-1,j))+1/dx^2;
        A(p(i+1,j),p(i+1,j))=A(p(i+1,j),p(i+1,j))+1/dx^2;
        % modify j-1,j+1 coefficients
        A(p(i,j-1),p(i,j-1))=A(p(i,j-1),p(i,j-1))+1/dy^2;
        A(p(i,j+1),p(i,j+1))=A(p(i,j+1),p(i,j+1))+1/dy^2;
    end
end

for i=nxb
    for j=nyb
        % zero out solid cells
        A(p(i,j),:)=0;
        A(p(i,j),p(i,j))=-2/dx^2-2/dy^2;
        b_vec(p(i,j))=0;
    end
end

% solve linear system

% % method 1: brute force
% phi_vec = A\b_vec;

% % method 2: LU
% [L,U] = lu(A);
% eta_vec = L\b_vec;
% phi_vec = U\eta_vec;

% % method 3: sorted LU
% [b_new,I] = sort(b_vec);
% [X,J] = sort(I); % J stores the inverse of sort such that b_vec = b_new(J)

% A_new = A(I,I);
% [L,U] = lu(A_new);
% eta_new = L\b_new;
% phi_new = U\eta_new;
% phi_vec = phi_new(J);

% method 4: sorted LU, minimal forward substitution
[b_new,I] = sort(b_vec);
[X,J] = sort(I); % J stores the inverse of sort such that b_vec = b_new(J)
m = find(b_new>0,1) 

n = nx*ny
A_new = A(I,I);
[L,U] = lu(A_new);
eta_new = zeros(n,1);
% minimal forward substitution, equivalent to eta_new = L\b_new
for j=m:n
    if L(j,j)==0
        disp(['stop: L matrix is singular at j=',num2str(j)])
        return
    end
    eta_new(j) = b_new(j)/L(j,j);
    iL=find( L((j+1):n,j)~=0 )+j;
    %for i=(j+1):n
    for i=iL
        b_new(i)=b_new(i)-L(i,j)*eta_new(j);
    end
end
% minimal backward substitution (expensive part), equivalent to phi_new = U\eta_new;
for j=n:-1:1
    if U(j,j)==0
        disp(['stop: U matrix is singular at j=',num2str(j)])
        return
    end
    phi_new(j) = eta_new(j)/U(j,j);
    iU=find( U(1:j-1,j)~=0 );
    %for i=1:(j-1)
    for i=iU
        eta_new(i) = eta_new(i)-U(i,j)*phi_new(j);
    end
end
phi_vec = phi_new(J);

% map phi to surf

for i=1:nx
    for j=1:ny
        phi(i,j) = phi_vec(p(i,j));
    end
end

figure
surf(x,y,phi')
xlabel('x')
ylabel('y')
zlabel('phi')






















