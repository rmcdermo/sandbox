% McDermott
% 10-29-2014
% build_sparse_matrix_block.m

function A = build_sparse_matrix_block(n,dx,bc)

% boundary condition type
% -----------------------
% bc(1) = left
% bc(2) = right
% bc(3) = top
% bc(4) = bottom
% values: 0 = Neumann, 1 = Dirichlet (ghost), 2 = Dirichlet (face)
%
% Consider that we are solving the following discrete equation in 2D:
%
% (p(i+1,j)-2*p(i,j)+p(i-1,j))/dx^2 + (p(i,j+1)-2*p(i,j)+p(i,j-1))^2 = b(i,j)
%
% Example of Neumann bc:
% Suppose we have dp/dx = Fx on the right-side bc.  The we have
% p(i+1)-p(i,j)= Fx*dx, and the eqn is rewritten as
%
% ( Fx*dx - p(i,j)+p(i-1,j) )/dx^2 + ... = b(i,j).
%
% In this case, the coefficient for p(i+1,j) is 0, and the coefficient
% for p(i,j) is -1.  And the source is augmented to be
% b(i,j) - Fx/dx.

nx = n(1);
ny = n(2);

dxdx = dx(1)^2;
dydy = dx(2)^2;

% Build A
A = sparse(nx*ny,nx*ny);
% interior nodes
for j = 2:ny-1
    for i = 2:nx-1
        p = (j-1)*nx+i;
        A(p,p+1) = 1/dxdx;
        A(p,p-1) = 1/dxdx;
        A(p,p+nx) = 1/dydy;
        A(p,p-nx) = 1/dydy;
        A(p,p) = -( 2/dxdx + 2/dydy );
    end
end
% left boundary
i = 1;
for j = 2:ny-1
    p = (j-1)*nx+i;
    A(p,p+1) = 1/dxdx;
    A(p,p+nx) = 1/dydy;
    A(p,p-nx) = 1/dydy;
    A(p,p) = -( (1+bc(1))/dxdx + 2/dydy );
end
% right boundary
i = nx;
for j = 2:ny-1
    p = (j-1)*nx+i;
    A(p,p-1) = 1/dxdx;
    A(p,p+nx) = 1/dydy;
    A(p,p-nx) = 1/dydy;
    A(p,p) = -( (1+bc(2))/dxdx + 2/dydy );
end
% top boundary
j = ny;
for i = 2:nx-1
    p = (j-1)*nx+i;
    A(p,p+1) = 1/dxdx;
    A(p,p-1) = 1/dxdx;
    A(p,p-nx) = 1/dydy;
    A(p,p) = -( 2/dxdx + (1+bc(3))/dydy );
end
% bottom boundary
j = 1;
for i = 2:nx-1
    p = (j-1)*nx+i;
    A(p,p+1) = 1/dxdx;
    A(p,p-1) = 1/dxdx;
    A(p,p+nx) = 1/dydy;
    A(p,p) = -( 2/dxdx + (1+bc(4))/dydy );
end

% left, bottom corner
i = 1; j = 1; p = (j-1)*nx+i;
A(p,p+1) = 1/dxdx;
A(p,p+nx) = 1/dydy;
A(p,p) = -( (1+bc(1))/dxdx + (1+bc(4))/dydy );

% left, top corner
i = 1; j = ny; p = (j-1)*nx+i;
A(p,p+1) = 1/dxdx;
A(p,p-nx) = 1/dydy;
A(p,p) = -( (1+bc(1))/dxdx + (1+bc(3))/dydy );

% right, bottom corner
i = nx; j = 1; p = (j-1)*nx+i;
A(p,p-1) = 1/dxdx;
A(p,p+nx) = 1/dydy;
A(p,p) = -( (1+bc(2))/dxdx + (1+bc(4))/dydy );

% right, top corner
i = nx; j = ny; p = (j-1)*nx+i;
A(p,p-1) = 1/dxdx;
A(p,p-nx) = 1/dydy;
A(p,p) = -( (1+bc(2))/dxdx + (1+bc(3))/dydy );

