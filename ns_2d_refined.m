% McDermott
% 8-24-15
% ns_2d_refined.m

close all
clear all

B = 2;       % solution amplitude
nu = 0;      % kinematic viscosity
Lx = 2*pi;   % x domain length, m
Ly = 2*pi;   % y domain length, m
T  = 2*pi;   % total time, s
cfl = 0.25;  % CFL

M(1).xs = 0;
M(1).xf = Lx;
M(1).ys = 0;
M(1).yf = Ly;
M(1).nx = 32;
M(1).ny = 32;
M(1).dx = (M(1).xf-M(1).xs)/M(1).nx;
M(1).dy = (M(1).yf-M(1).ys)/M(1).ny;

M(2).xs = Lx/4;
M(2).xf = Lx-Lx/4;
M(2).ys = Ly/4;
M(2).yf = Ly-Ly/4;
M(2).nx = 32;
M(2).ny = 32;
M(2).dx = (M(2).xf-M(2).xs)/M(2).nx;
M(2).dy = (M(2).yf-M(2).ys)/M(2).ny;

M(1).x = M(1).xs + [0:M(1).nx]*M(1).dx;
M(1).y = M(1).ys + [0:M(1).ny]*M(1).dy;
M(2).x = M(2).xs + [0:M(2).nx]*M(2).dx;
M(2).y = M(2).ys + [0:M(2).ny]*M(2).dy;

M(1).xc = M(1).x(1:M(1).nx) + 0.5*M(1).dx;
M(1).yc = M(1).y(1:M(1).ny) + 0.5*M(1).dy;

M(2).xc = M(2).x(1:M(2).nx) + 0.5*M(2).dx;
M(2).yc = M(2).y(1:M(2).ny) + 0.5*M(2).dy;

% generate initial conditions (staggered, face-low storage)

% for a given cell...

%                     ^ v(i,j+1)
%            ---------|--------
%            |                |
%            |                |
%            |                |
%    u(i,j) --->   H(i,j)    ---> u(i+1,j)
%            |                |
%            |                |
%            |       ^        |
%            --------|---------
%                      v(i,j)

t = 0;
dt = cfl*M(2).dx/(B+1);
for nm = 1:2
    for j=1:M(nm).ny
        for i=1:M(nm).nx
            M(nm).p(i,j) = -(B^2)/4*( cos(2*(M(nm).xc(i)-t)) + cos(2*(M(nm).yc(j)-t)) )*exp(-4*nu*t);
        end
    end
    for j=1:M(nm).ny
        for i=1:M(nm).nx+1
            M(nm).u(i,j) = 1 - B*cos(M(nm).x(i)-t)*sin(M(nm).yc(j)-t)*exp(-2*nu*t);
        end
    end
    for j=1:M(nm).ny+1
        for i=1:M(nm).nx
            M(nm).v(i,j) = 1 + B*sin(M(nm).xc(i)-t)*cos(M(nm).y(j)-t)*exp(-2*nu*t);
        end
    end
end

% % plot rectangular array of cells
% pcolor(M(1).x,M(1).y,pad(M(1).p)); hold on
% pcolor(M(2).x,M(2).y,pad(M(2).p))
% axis square
% %colormap(hot(20))

% % compute finite volume (integral) cell divergence
% for nm = 1:2
%     for j=1:M(nm).ny
%         for i=1:M(nm).nx
%             M(nm).D(i,j) = ( M(nm).u(i+1,j)-M(nm).u(i,j) ) * M(nm).dy + ( M(nm).v(i,j+1)-M(nm).v(i,j) ) * M(nm).dx;
%         end
%     end
% end

% % plot cell divergence
% figure
% pcolor(M(1).x,M(1).y,pad(M(1).D)); hold on
% pcolor(M(2).x,M(2).y,pad(M(2).D))
% axis square
% return

% see locate_mesh.m for description of i_lo, i_hi, etc.
[index_list,map1,map2,imap,jmap,cells,rx,ry,n_cells,ierror] = locate_mesh(M(1),M(2)); 
i_lo = index_list(1);
i_hi = index_list(2);
j_lo = index_list(3);
j_hi = index_list(4);
%return

% build A matrix

A = sparse(n_cells,n_cells)

for row=1:n_cells

    % east face coefficients

    DX = 0.5 * ( M(cells(row).east_level).dx + M(cells(row).level).dx );
    DY = min( M(cells(row).level).dy, M(cells(row).east_level).dy );
    NJ = find(cells(row).east_index>0);
    for jj=NJ
        A(row,row) = A(row,row) - 1/DX * DY;
        A(row,cells(row).east_index(jj)) = 1/DX * DY;
    end

    % west face coefficients

    DX = 0.5 * ( M(cells(row).west_level).dx + M(cells(row).level).dx );
    DY = min( M(cells(row).level).dy, M(cells(row).west_level).dy );
    NJ = find(cells(row).west_index>0);
    for jj=NJ
        A(row,row) = A(row,row) - 1/DX * DY;
        A(row,cells(row).west_index(jj)) = 1/DX * DY;
    end

    % north face coefficients

    DY = 0.5 * ( M(cells(row).north_level).dy + M(cells(row).level).dy );
    DX = min( M(cells(row).level).dx, M(cells(row).north_level).dx );
    NI = find(cells(row).north_index>0);
    for ii=NI
        A(row,row) = A(row,row) - 1/DY * DX;
        A(row,cells(row).north_index(ii)) = 1/DY * DX;
    end

    % south face coefficients

    DY = 0.5 * ( M(cells(row).south_level).dy + M(cells(row).level).dy );
    DX = min( M(cells(row).level).dx, M(cells(row).south_level).dx );
    NI = find(cells(row).south_index>0);
    for ii=NI
        A(row,row) = A(row,row) - 1/DY * DX;
        A(row,cells(row).south_index(ii)) = 1/DY * DX;
    end

end

% % uncomment to view matrix structure

% figure
% full(A)
% spy(A)
% eig(full(A))
% tightfig(gcf);
% print(gcf,'-dpdf','spyA')
% return

while t<T % time loop
    
    t = t + dt;

    % do coarse grid stresses (periodic bcs)
    for j=1:M(1).ny
        for i=1:M(1).nx
            
            ip1=i+1;
            im1=i-1;
            jp1=j+1;
            jm1=j-1;
        
            if ip1>M(1).nx; ip1=ip1-M(1).nx; end
            if jp1>M(1).ny; jp1=jp1-M(1).ny; end
            if im1<1;  im1=im1+M(1).nx; end
            if jm1<1;  jm1=jm1+M(1).ny; end
             
            % compute advective stress (Stokes form)
            M(1).omega3(i,j) = (M(1).v(i,j)-M(1).v(im1,j))/M(1).dx - (M(1).u(i,j)-M(1).u(i,jm1))/M(1).dy;
            M(1).ubar(i,j) = 0.5*(M(1).u(i,jm1)+M(1).u(i,j));
            M(1).vbar(i,j) = 0.5*(M(1).v(im1,j)+M(1).v(i,j));
            
            % compute viscous stress components
            dudy = (M(1).u(i,j)-M(1).u(i,jm1))/M(1).dy;
            dvdx = (M(1).v(i,j)-M(1).v(im1,j))/M(1).dx;
            dudx = (M(1).u(ip1,j)-M(1).u(i,j))/M(1).dx;
            dvdy = (M(1).v(i,jp1)-M(1).v(i,j))/M(1).dy;

            M(1).tau11(i,j) = -2*nu*dudx; % cell center
            M(1).tau22(i,j) = -2*nu*dvdy; % cell center
            M(1).tau12(i,j) = -nu*(dudy + dvdx); % vertex

        end
    end
    
    % compute force terms and update velocity predictor on coarse grid
    for j=1:M(1).ny
        for i=1:M(1).nx
            
            ip1=i+1;
            im1=i-1;
            jp1=j+1;
            jm1=j-1;
            
            if ip1>M(1).nx; ip1=ip1-M(1).nx; end
            if jp1>M(1).ny; jp1=jp1-M(1).ny; end
            if im1<1;  im1=im1+M(1).nx; end
            if jm1<1;  jm1=jm1+M(1).ny; end

            % Stokes form
            Fx = -0.5*( M(1).vbar(i,j)*M(1).omega3(i,j) + M(1).vbar(i,jp1)*M(1).omega3(i,jp1) );
            Fy =  0.5*( M(1).ubar(i,j)*M(1).omega3(i,j) + M(1).ubar(ip1,j)*M(1).omega3(ip1,j) );

            % add viscous stresses
            Fx = Fx + (M(1).tau11(i,j)-M(1).tau11(im1,j))/M(1).dx + (M(1).tau12(ip1,j)-M(1).tau12(i,j))/M(1).dy;
            Fy = Fy + (M(1).tau12(ip1,j)-M(1).tau12(i,j))/M(1).dx + (M(1).tau22(i,j)-M(1).tau22(i,jm1))/M(1).dy;
            
            % velocity predictor
            M(1).uhat(i,j) = M(1).u(i,j) - dt*Fx;
            M(1).vhat(i,j) = M(1).v(i,j) - dt*Fy;
        end
    end

    % apply periodic bcs to uhat, vhat on coarse grid
    M(1).uhat(M(1).nx+1,:) = M(1).uhat(1,:);
    M(1).vhat(:,M(1).ny+1) = M(1).vhat(:,1);

    % do fine grid diagonal stress components (embedded bcs), stored at cell center
    for jj=1:M(2).ny
        for ii=1:M(2).nx

            UP = M(2).u(ii+1,jj);
            UM = M(2).u(ii,jj);
            DX = M(2).dx;

            VP = M(2).v(ii,jj+1);
            VM = M(2).v(ii,jj);
            DY = M(2).dy;

            % compute viscous stress components
            M(2).tau11(ii,jj) = -2*nu*(UP-UM)/DX; % cell center
            M(2).tau22(ii,jj) = -2*nu*(VP-VM)/DY; % cell center
        end
    end

    % do fine grid off-diagonal stress components (embedded bcs), stored at vertex, cell-low
    for jj=1:M(2).ny+1
        for ii=1:M(2).nx+1

            DX = M(2).dx;

            if ii<=M(2).nx
                VXP = M(2).v(ii,jj);
            else
                % linear interpolation for guard cell velocity from coarse mesh
                % this probably requires some explanation
                j = j_lo - 1 + ceil(jj/ry);
                Y = [M(1).y(j) M(1).y(j+1)];
                V = [M(1).v(i_hi+1,j) M(1).v(i_hi+1,j+1)];
                VXP = interp1(Y,V,M(2).y(jj));
                DX = 0.5 * (M(1).dx + M(2).dx);
            end

            if ii>1
                VXM = M(2).v(ii-1,jj);
            else
                j = j_lo - 1 + ceil(jj/ry);
                Y = [M(1).y(j) M(1).y(j+1)];
                V = [M(1).v(i_lo-1,j) M(1).v(i_lo-1,j+1)];
                VXM = interp1(Y,V,M(2).y(jj));
                DX = 0.5 * (M(1).dx + M(2).dx);
            end

            DY = M(2).dy;

            if jj<=M(2).ny
                UYP = M(2).u(ii,jj);
            else
                i = i_lo - 1 + ceil(ii/ry);
                X = [M(1).x(i) M(1).x(i+1)];
                U = [M(1).u(i,j_hi+1) M(1).u(i+1,j_hi+1)];
                UYP = interp1(X,U,M(2).x(ii));
                DY = 0.5 * (M(1).dy + M(2).dy);
            end

            if jj>1
                UYM = M(2).u(ii,jj-1);
            else
                i = i_lo - 1 + ceil(ii/ry);
                X = [M(1).x(i) M(1).x(i+1)];
                U = [M(1).u(i,j_lo-1) M(1).u(i+1,j_lo-1)];
                UYM = interp1(X,U,M(2).x(ii));
                DY = 0.5 * (M(1).dy + M(2).dy);
            end

            dudy = (UYP-UYM)/DY;
            dvdx = (VXP-VXM)/DX;

            % compute advective stress (Stokes form)
            M(2).omega3(ii,jj) = dvdx - dudy;
            M(2).ubar(ii,jj) = 0.5*(UYM+UYP);
            M(2).vbar(ii,jj) = 0.5*(VXM+VXP);
            
            % compute viscous stress components
            M(2).tau12(ii,jj) = -nu*(dudy + dvdx); % vertex
        end
    end

    % compute force terms and update velocity predictor on fine grid -- x
    for jj=1:M(2).ny
        for ii=1:M(2).nx+1
            % Stokes form
            Fx = -0.5*( M(2).vbar(ii,jj)*M(2).omega3(ii,jj) + M(2).vbar(ii,jj+1)*M(2).omega3(ii,jj+1) );

            % add viscous stresses
            DX = M(2).dx;
            if ii<=M(2).nx
                T11_P = M(2).tau11(ii,jj);
            else
                j = j_lo - 1 + ceil(jj/ry);
                T11_P = M(1).tau11(i_hi+1,j);
            end
            if ii>1
                T11_M = M(2).tau11(ii-1,jj);
            else
                j = j_lo - 1 + ceil(jj/ry);
                T11_M = M(1).tau11(i_lo-1,j);
            end

            Fx = Fx + (T11_P-T11_M)/DX + (M(2).tau12(ii,jj+1)-M(2).tau12(ii,jj))/M(2).dy;
            
            % velocity predictor
            M(2).uhat(ii,jj) = M(2).u(ii,jj) - dt*Fx;
        end
    end

    % compute force terms and update velocity predictor on fine grid -- y

    for jj=1:M(2).ny+1
        for ii=1:M(2).nx
            % Stokes form
            Fy =  0.5*( M(2).ubar(ii,jj)*M(2).omega3(ii,jj) + M(2).ubar(ii+1,jj)*M(2).omega3(ii+1,jj) );

            % add viscous stresses
            DY = M(2).dy;
            if jj<=M(2).ny
                T22_P = M(2).tau22(ii,jj);
            else
                i = i_lo - 1 + ceil(ii/rx);
                T22_P = M(1).tau22(i,j_hi+1);
            end

            if jj>1
                T22_M = M(2).tau22(ii,jj-1);
            else
                i = i_lo - 1 + ceil(ii/rx);
                T22_M = M(1).tau22(i,j_lo-1);
            end

            Fy = Fy + (M(2).tau12(ii+1,jj)-M(2).tau12(ii,jj))/M(2).dx + (T22_P-T22_M)/DY;
            
            % velocity predictor
            M(2).vhat(ii,jj) = M(2).v(ii,jj) - dt*Fy;
        end
    end

    % correct cells on coarse-fine interface

    % this loop does all coarse cells under refinement region and along coarse-fine boundary
    for i=i_lo:i_hi+1
        ii = (i-i_lo)*rx + 1;
        for j=j_lo:j_hi
            jj = (j-j_lo)*ry*ones(1,ry) + [1:ry];
            M(1).uhat(i,j) = mean( M(2).uhat(ii,jj) );
        end
    end

    for j = j_lo:j_hi+1
        jj = (j-j_lo)*ry + 1;
        for i=i_lo:i_hi
            ii = (i-i_lo)*rx*ones(1,rx) + [1:rx];
            M(1).vhat(i,j) = mean( M(2).vhat(ii,jj) );
        end
    end

    % build right-hand-side of Poisson equation

    for nm=1:2
        for j=1:M(nm).nx
            for i=1:M(nm).ny
                M(nm).D(i,j) = ( M(nm).uhat(i+1,j) - M(nm).uhat(i,j) ) * M(nm).dy ...
                             + ( M(nm).vhat(i,j+1) - M(nm).vhat(i,j) ) * M(nm).dx;
            end
        end
    end

    % assemble cell divergence into RHS vector

    for row=1:n_cells
        b(row) = M(cells(row).level).D(imap(row),jmap(row));
    end

    % solve Poisson equation for "condensed" system
    Ac = A(1:n_cells-1,1:n_cells-1);
    bc = b(1:n_cells-1);
    pc = Ac\bc';
    pvec = [pc;0];
    
    % map solution vector to computational indices
    for row=1:n_cells
        M(cells(row).level).p(imap(row),jmap(row)) = pvec(row);
    end

    % project velocities
    % note: dt may be omitted here if it is left out of the b vector (right hand side) above

    % coarse grid

    for j=1:M(1).ny
        for i=1:M(1).nx
            im1=i-1;
            jm1=j-1;
            if im1<1; im1=im1+M(1).nx; end
            if jm1<1; jm1=jm1+M(1).ny; end
            
            M(1).u(i,j) = M(1).uhat(i,j) - ( M(1).p(i,j) - M(1).p(im1,j) ) / M(1).dx;
            M(1).v(i,j) = M(1).vhat(i,j) - ( M(1).p(i,j) - M(1).p(i,jm1) ) / M(1).dy;
        end
    end

    % apply periodic bcs to u, v on coarse grid
    M(1).u(M(1).nx+1,:) = M(1).u(1,:);
    M(1).v(:,M(1).ny+1) = M(1).v(:,1);

    % fine grid, interior unknowns

    for jj=1:M(2).ny
        for ii=2:M(2).nx
            M(2).u(ii,jj) = M(2).uhat(ii,jj) - ( M(2).p(ii,jj) - M(2).p(ii-1,jj) ) / M(2).dx;
        end
    end

    for jj=2:M(2).ny
        for ii=1:M(2).nx
            M(2).v(ii,jj) = M(2).vhat(ii,jj) - ( M(2).p(ii,jj) - M(2).p(ii,jj-1) ) / M(2).dy;
        end
    end

    % fine grid east side

    ii = M(2).nx+1;
    DX = 0.5 * ( M(1).dx + M(2).dx );
    for jj=1:M(2).ny
        j = j_lo - 1 + ceil(jj/ry);
        M(2).u(ii,jj) = M(2).uhat(ii,jj) - ( M(1).p(i_hi+1,j) - M(2).p(ii-1,jj) ) / DX;
    end

    % fine grid west side

    ii = 1;
    DX = 0.5 * ( M(1).dx + M(2).dx );
    for jj=1:M(2).ny
        j = j_lo - 1 + ceil(jj/ry);
        M(2).u(ii,jj) = M(2).uhat(ii,jj) - ( M(2).p(ii,jj) - M(1).p(i_lo-1,j) ) / DX;
    end

    % fine grid north side

    jj = M(2).ny + 1;
    DY = 0.5 * ( M(1).dy + M(2).dy );
    for ii=1:M(2).nx
        i = i_lo - 1 + ceil(ii/rx);
        M(2).v(ii,jj) = M(2).vhat(ii,jj) - ( M(1).p(i,j_hi+1) - M(2).p(ii,jj-1) ) / DY;
    end

    % fine grid south side

    jj = 1;
    DY = 0.5 * ( M(1).dy + M(2).dy );
    for ii=1:M(2).nx
        i = i_lo - 1 + ceil(ii/rx);
        M(2).v(ii,jj) = M(2).vhat(ii,jj) - ( M(2).p(ii,jj) - M(1).p(i,j_lo-1) ) / DY;
    end

    % correct velocity components on coarse grid, restrict fine grid components

    % this loop does all coarse cells under refinement region and along coarse-fine boundary
    for i=i_lo:i_hi+1
        ii = (i-i_lo)*rx + 1;
        for j=j_lo:j_hi
            jj = (j-j_lo)*ry*ones(1,ry) + [1:ry];
            M(1).u(i,j) = mean( M(2).u(ii,jj) );
        end
    end

    for j = j_lo:j_hi+1
        jj = (j-j_lo)*ry + 1;
        for i=i_lo:i_hi
            ii = (i-i_lo)*rx*ones(1,rx) + [1:rx];
            M(1).v(i,j) = mean( M(2).v(ii,jj) );
        end
    end

    % check divergence

    for nm=1:2
        for j=1:M(nm).nx
            for i=1:M(nm).ny
                M(nm).D(i,j) = ( M(nm).u(i+1,j) - M(nm).u(i,j) ) * M(nm).dy ...
                             + ( M(nm).v(i,j+1) - M(nm).v(i,j) ) * M(nm).dx;
            end
        end
    end

    % % plot divergence
    % hold off
    % pcolor(M(1).x,M(1).y,pad(M(1).D)); hold on
    % pcolor(M(2).x,M(2).y,pad(M(2).D))
    % axis square

    display( ['max div M1 = ',num2str(max(max(abs(M(1).D)))),'   max div M2 = ',num2str(max(max(abs(M(2).D))))] )

    subplot(1,2,1), surf(M(1).x,M(1).yc,M(1).u'), axis([0 Lx 0 Ly -1 3])
    subplot(1,2,2), surf(M(2).x,M(2).yc,M(2).u'), axis([0 Lx 0 Ly -1 3])
    
    pause(0.001)
    %return

end % time loop




















