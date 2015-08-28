% McDermott
% 10-29-2014
% Laplace_2d_block.m
%
% Consider the domain D formed by combining a set of blocks such that an internal boundary exists.
% We want to generate a matrix for the Laplace equation on D = sum(Di) that excludes the immersed
% boundary block B.
%
%                -------------------------------
%                |         :         :         |
%                |    D6   :    D7   :    D8   |
%                |         :         :         |
%                |         :         :         |
%                | - - - - ----------- - - - - |
%                |         |BBBBBBBBB|         |
%                |    D4   |BBBBBBBBB|    D5   |
%                |         |BBBBBBBBB|         |
%                |         |BBBBBBBBB|         |
%                | - - - - ----------- - - - - |
%                |         :         :         |
%                |    D1   :    D2   :    D3   |
%                |         :         :         |
%                |         :         :         |
%                -------------------------------

close all
clear all

nb = 8; % number of blocks

xs = [0 1 2 0 2 0 1 2]; % x start for each block
xf = [1 2 3 1 3 1 2 3]; % x finish for each block
ys = [0 0 0 1 1 2 2 2]; % y start for each block
yf = [1 1 1 2 2 3 3 3]; % y finish for each block
nx = [8 8 8 8 8 8 8 8];
ny = nx;

for ib=1:nb
    D(ib).xf = xf(ib);
    D(ib).xs = xs(ib);
    D(ib).yf = yf(ib);
    D(ib).ys = ys(ib);
    D(ib).nx = nx(ib);
    D(ib).ny = ny(ib);
    D(ib).Lx = D(ib).xf - D(ib).xs;
    D(ib).Ly = D(ib).yf - D(ib).ys;
    D(ib).dx = D(ib).Lx/D(ib).nx;
    D(ib).dy = D(ib).Ly/D(ib).ny;
    D(ib).x = D(ib).xs+0.5*D(ib).dx:D(ib).dx:D(ib).xf-0.5*D(ib).dx;
    D(ib).y = D(ib).ys+0.5*D(ib).dy:D(ib).dy:D(ib).yf-0.5*D(ib).dy;
end

% define block neighbors

D(1).East_Neighbor = 2;
D(1).West_Neighbor = 0;
D(1).North_Neighbor = 4;
D(1).South_Neighbor = 0;

D(2).East_Neighbor = 3;
D(2).West_Neighbor = 1;
D(2).North_Neighbor = 0;
D(2).South_Neighbor = 0;

D(3).East_Neighbor = 0;
D(3).West_Neighbor = 2;
D(3).North_Neighbor = 5;
D(3).South_Neighbor = 0;

D(4).East_Neighbor = 0;
D(4).West_Neighbor = 0;
D(4).North_Neighbor = 6;
D(4).South_Neighbor = 1;

D(5).East_Neighbor = 0;
D(5).West_Neighbor = 0;
D(5).North_Neighbor = 8;
D(5).South_Neighbor = 3;

D(6).East_Neighbor = 7;
D(6).West_Neighbor = 0;
D(6).North_Neighbor = 0;
D(6).South_Neighbor = 4;

D(7).East_Neighbor = 8;
D(7).West_Neighbor = 6;
D(7).North_Neighbor = 0;
D(7).South_Neighbor = 0;

D(8).East_Neighbor = 0;
D(8).West_Neighbor = 7;
D(8).North_Neighbor = 0;
D(8).South_Neighbor = 5;

% build each block
%
% boundary condition type
% -----------------------
% bc(1) = left
% bc(2) = right
% bc(3) = top
% bc(4) = bottom
% values: 0 = Neumann, 1 = Dirichlet (ghost), 2 = Dirichlet (face)
D(1).A = build_sparse_matrix_block([D(1).nx D(1).ny],[D(1).dx D(1).dy],[2 1 1 0]);
D(2).A = build_sparse_matrix_block([D(2).nx D(2).ny],[D(2).dx D(2).dy],[1 1 0 0]);
D(3).A = build_sparse_matrix_block([D(3).nx D(3).ny],[D(3).dx D(3).dy],[1 2 1 0]);
D(4).A = build_sparse_matrix_block([D(4).nx D(4).ny],[D(4).dx D(4).dy],[2 0 1 1]);
D(5).A = build_sparse_matrix_block([D(5).nx D(5).ny],[D(5).dx D(5).dy],[0 2 1 1]);
D(6).A = build_sparse_matrix_block([D(6).nx D(6).ny],[D(6).dx D(6).dy],[2 1 0 1]);
D(7).A = build_sparse_matrix_block([D(7).nx D(7).ny],[D(7).dx D(7).dy],[1 1 0 0]);
D(8).A = build_sparse_matrix_block([D(8).nx D(8).ny],[D(8).dx D(8).dy],[1 2 0 1]);

% assemble A
A = sparse([]);
for ib=1:nb
    S = size(A);
    D(ib).P0 = S(1);
    N = D(ib).nx*D(ib).ny;
    A(S(1)+1:S(1)+N,S(2)+1:S(2)+N) = D(ib).A;
    % create mapping for ib,i,j
    for i=1:D(ib).nx
        for j=1:D(ib).ny
            P(ib,i,j) = D(ib).P0 + (j-1)*D(ib).nx + i;
        end
    end
end

% fill in block neighbor elements
for ib=1:nb
    % fill in East Neighbor
    ibn=D(ib).East_Neighbor;
    if ibn>0
        for j=1:D(ib).ny
            A(P(ib,D(ib).nx,j),P(ibn,1,j)) = 1/D(ib).dx^2;
        end
    end
    % fill in West Neighbor
    ibn=D(ib).West_Neighbor;
    if ibn>0
        for j=1:D(ib).ny
            A(P(ib,1,j),P(ibn,D(ibn).nx,j)) = 1/D(ib).dx^2;
        end
    end
    % fill in North Neighbor
    ibn=D(ib).North_Neighbor;
    if ibn>0
        for i=1:D(ib).nx
            A(P(ib,i,D(ib).ny),P(ibn,i,1)) = 1/D(ib).dy^2;
        end
    end
    % fill in South Neighbor
    ibn=D(ib).South_Neighbor;
    if ibn>0
        for i=1:D(ib).nx
            A(P(ib,i,1),P(ibn,i,D(ibn).ny)) = 1/D(ib).dy^2;
        end
    end
end

% source term for each block
for ib=1:nb
    D(ib).b = zeros(D(ib).nx*D(ib).ny,1);
end

% left boundary blocks
for ib=[1,4,6]
    i=1; % left boundary of block ib
    for j = 1:D(ib).ny
        p = (j-1)*D(ib).nx+i;
        D(ib).b(p) = -2/D(ib).dx^2;
    end
end

% assemble b
b = sparse([]);
for ib=1:nb
    N = D(ib).nx*D(ib).ny;
    b(D(ib).P0+1:D(ib).P0+N) = D(ib).b;
end

% LU decomposition of A

[L,U] = lu(A);

eta_vec = L\b';
phi_vec = U\eta_vec;

% map phi to surf

for ib=1:nb
    for i=1:D(ib).nx
        for j=1:D(ib).ny
            D(ib).phi(i,j) = phi_vec(P(ib,i,j));
        end
    end
    surf(D(ib).x,D(ib).y,D(ib).phi'); hold on
    %pause
end

% surf(x,y,phi')
% xlabel('x')
% ylabel('y')
% zlabel('\phi')





