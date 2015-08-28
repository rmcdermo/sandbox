% McDermott
% 12-3-2014
% build_sparse_matrix_2.m

function A = build_sparse_matrix_2(n,dx,bc)

% boundary condition type
% -----------------------
% bc(1) = west
% bc(2) = east
% bc(3) = south
% bc(4) = north
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
for j = 1:ny
    for i = 1:nx
        p = (j-1)*nx+i;

        % west
        if i>1
            bcw = 1;
            A(p,p-1) = bcw/dxdx;
        else
            bcw = bc(1);
        end

        % east
        if i<nx
            bce = 1;
            A(p,p+1) = bce/dxdx;
        else
            bce = bc(2);
        end

        % south
        if j>1
            bcs = 1;
            A(p,p-nx) = bcs/dydy;
        else
            bcs = bc(3);
        end

        % north
        if j<ny
            bcn = 1;
            A(p,p+nx) = bcn/dydy;
        else
            bcn = bc(4);
        end
        
        A(p,p) = -( (bcw+bce)/dxdx + (bcs+bcn)/dydy );
    end
end

