% McDermott
% 8-25-15
% locate_mesh.m
% 
% Translated from SUBROUTINE LOCATE_MESH in samr.f90 in FDS.
%
%   Assumes uniform grid in each direction.
%
%   -------------------------------
%   |         |         |         |
%   |         |         |         |<---MESHES(M1)
%   |         |         |         |
%   |         |         |         |
%   -------------------------------
%   |         |-|-|-|-|-|         |
%   |         |-|-|-|-|-|<-------------MESHES(M2)
%   |         |-|-|-|-|-|         |
%   |         |-|-|-|-|-|         |
%   -------------------------------
%   |         |         |         |
%   |         |         |         |
%   |         |         |         |
%   |         |         |         |
%   -------------------------------

function [INDEX_LIST,MAP1,MAP2,IMAP,JMAP,C,rx,ry,n_cells,IERROR] = locate_mesh(M1,M2)

IERROR=0;
INDEX_LIST=0;

% M1 => coarse mesh
% M2 => fine mesh

% Locate fine mesh within coarse mesh.
% In coarse mesh indices, the fine mesh spans I_LO to I_HI in x and J_LO to J_HI in y.

I_LO = max(1,     round((M2.xs-M1.xs)/M1.dx)+1 );
I_HI = min(M1.nx, round((M2.xf-M1.xs)/M1.dx)   );
if I_LO>M1.nx | I_HI<1 % meshes do not overlap
  IERROR=1;
  return
end

J_LO = max(1,     round((M2.ys-M1.ys)/M1.dy)+1 );
J_HI = min(M1.ny, round((M2.yf-M1.ys)/M1.dy)   );
if J_LO>M1.ny | J_HI<1 % meshes do not overlap
  IERROR=1;
  return
end

% Find fine mesh off-set.
% If the fine mesh is not completely contained within the coarse mesh, the first fine mesh
% index on the coarse mesh is II_LO in x and JJ_LO in y.

II_LO = max(0, round((M1.xs-M2.xs)/M2.dx) );
JJ_LO = max(0, round((M1.ys-M2.ys)/M2.dy) );

INDEX_LIST = [I_LO,I_HI,J_LO,J_HI,II_LO,JJ_LO];

% Map the mesh indices to the global solution vector

counter=1; % index of global solution vector

% do coarse mesh map, MAP1
for j=1:M1.ny
    for i=1:M1.nx
        if i>=I_LO & i<=I_HI & j>=J_LO & j<=J_HI
            MAP1(i,j)=0;
        else
            MAP1(i,j)=counter;
            IMAP(counter)=i;
            JMAP(counter)=j;
            counter=counter+1;
        end
    end
end

n_coarse_cells = counter-1;

% do fine mesh map, MAP2
for jj=1:M2.ny
    for ii=1:M2.nx
        MAP2(ii,jj)=counter;
        IMAP(counter)=ii;
        JMAP(counter)=jj;
        counter=counter+1;
    end
end

n_cells = counter-1;

% refinement ratio

rx = M1.dx/M2.dx;
ry = M1.dy/M2.dy;

% initialize neighbor slots

C.east_index = zeros(1,ry);
C.west_index = zeros(1,ry);
C.north_index = zeros(1,rx);
C.south_index = zeros(1,rx);

% Find global index of neighbors for each cell

for n = 1:n_cells % n_cells_for

    if n<=n_coarse_cells
        C(n).level=1;
    else
        C(n).level=2;
    end

    if C(n).level==1 % level_1_if
        i = IMAP(n);
        j = JMAP(n);

        % for coarse grid we need to handle periodic bcs
        ip1 = i+1;
        im1 = i-1;
        jp1 = j+1;
        jm1 = j-1;
        if ip1>M1.nx; ip1=ip1-M1.nx; end
        if im1<1    ; im1=im1+M1.nx; end
        if jp1>M1.ny; jp1=jp1-M1.ny; end
        if jm1<1    ; jm1=jm1+M1.ny; end

        % east neighbor(s)
        if MAP1(ip1,j)==0
            C(n).east_level = 2;
            for nn=1:ry
                jj = (j-J_LO)*ry + nn;
                C(n).east_index(nn) = MAP2(1,jj);
            end
        else
            C(n).east_level = 1;
            C(n).east_index(1) = MAP1(ip1,j);
        end

        % west neighbor(s)
        if MAP1(im1,j)==0
            C(n).west_level = 2;
            for nn=1:ry
                jj = (j-J_LO)*ry + nn;
                C(n).west_index(nn) = MAP2(M2.nx,jj);
            end
        else
            C(n).west_level = 1;
            C(n).west_index(1) = MAP1(im1,j);
        end

        % north neighbor(s)
        if MAP1(i,jp1)==0
            C(n).north_level = 2;
            for nn=1:rx
                ii = (i-I_LO)*rx + nn;
                C(n).north_index(nn) = MAP2(ii,1);
            end
        else
            C(n).north_level = 1;
            C(n).north_index(1) = MAP1(i,jp1);
        end

        % south neighbor(s)
        if MAP1(i,jm1)==0
            C(n).south_level = 2;
            for nn=1:rx
                ii = (i-I_LO)*rx + nn;
                C(n).south_index(nn) = MAP2(ii,M2.ny);
            end
        else
            C(n).south_level = 1;
            C(n).south_index(1) = MAP1(i,jm1);
        end

    elseif C(n).level==2

        ii = IMAP(n);
        jj = JMAP(n);

        % east neighbor
        if ii+1>M2.nx
            j = J_LO - 1 + ceil(jj/ry);
            C(n).east_level = 1;
            C(n).east_index(1) = MAP1(I_HI+1,j);
        else
            C(n).east_level = 2;
            C(n).east_index(1) = MAP2(ii+1,jj);
        end

        % west neighbor
        if ii-1<1
            j = J_LO - 1 + ceil(jj/ry);
            C(n).west_level = 1;
            C(n).west_index(1) = MAP1(I_LO-1,j);
        else
            C(n).west_level = 2;
            C(n).west_index(1) = MAP2(ii-1,jj);
        end

        % north neighbor
        if jj+1>M2.ny
            i = I_LO - 1 + ceil(ii/rx);
            C(n).north_level = 1;
            C(n).north_index(1) = MAP1(i,J_HI+1);
        else
            C(n).north_level = 2;
            C(n).north_index(1) = MAP2(ii,jj+1);
        end

        % south neighbor
        if jj-1<1
            i = I_LO - 1 + ceil(ii/rx);
            C(n).south_level = 1;
            C(n).south_index(1) = MAP1(i,J_LO-1);
        else
            C(n).south_level = 2;
            C(n).south_index(1) = MAP2(ii,jj-1);
        end

    end % level_1_if

end % n_cells_for

return





































