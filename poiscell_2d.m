% McDermott
% 10-14-2015
% poiscell_2d.m
%
% Solve Poisson equation on a 5-point Cartesian stencil with immersed cutcells.
%
%            ......---^---......
%            :     |     |     :
%            :     >  o  >     :
%            :     |     |     :
%            |--^--|-/^\-|--^--|
%            |     /     \     |
%            >  o  >  o  >  o  >
%            |     \     /     |
%            |--^--|-\^/-|--^--| 
%            :     |     |     :
%            :     >  o  >     :
%            :     |     |     :
%            ......---^---......
%
% Imagine that at the Cartesian level the face velocities obey the Cartesian divergence constraint.
% The scalar bcs around the outer boundary are homogeneous Neumann.  The goal is to reconstruct the
% velocities on the cutfaces such that the divergence constraint for all cutcell volumes are obeyed.
% If there is only one gas phase cutface, there is no issue with reconstruction---a simple rescaling
% of the face velocity component will do.  But here we consider that there may be two gas phase
% cutfaces to reconstruct.  Each cutcell volume stores a scalar pressure variable at its centroid
% for use in a projection.

close all
clear all

% create a 3 x 3 array of cells with a nontrivial solenoidal velocity field

poiscell_2d_init

umax = max(max(abs(u)));
vmax = max(max(abs(v)));

% draw mesh
vscal = 0.25; % velocity scale factor
i = [2 1 2 3 2];
j = [1 2 2 2 3];
for c=1:5
    x_vert = [x(i(c)) x(i(c)+1) x(i(c)+1) x(i(c))];
    y_vert = [y(j(c)) y(j(c)) y(j(c)+1) y(j(c)+1)];
    patch(x_vert,y_vert,'w'); hold on

    % while we're here, draw velocity arrows
    plot([x(i(c))-u(i(c),j(c))/umax*dx*vscal,x(i(c))+u(i(c),j(c))/umax*dx*vscal],[yp(j(c)),yp(j(c))],'r-');
    plot([x(i(c)+1)-u(i(c)+1,j(c))/umax*dx*vscal,x(i(c)+1)+u(i(c)+1,j(c))/umax*dx*vscal],[yp(j(c)),yp(j(c))],'r-');

    plot([xp(i(c)),xp(i(c))],[y(j(c))-v(i(c),j(c))/vmax*dy*vscal,y(j(c))+v(i(c),j(c))/vmax*dy*vscal],'b-');
    plot([xp(i(c)),xp(i(c))],[y(j(c)+1)-v(i(c),j(c)+1)/vmax*dy*vscal,y(j(c)+1)+v(i(c),j(c)+1)/vmax*dy*vscal],'b-');

end

% draw immersed boundary

x_vert = [x(2)+dx/2 x(3)+dx/4 x(2)+dx/2 x(2)-dx/4];
y_vert = [y(3)+dy/4 y(2)+dy/2 y(2)-dy/4 y(2)+dy/2];
patch(x_vert,y_vert,'b');

% interpolate vel components to cutfaces with zero at surface

us(1) = u(2,1)/6;
us(2) = u(2,3)/6;
us(3) = u(3,3)/6;
us(4) = u(3,1)/6;

vs(1) = v(1,2)/6;
vs(2) = v(1,3)/6;
vs(3) = v(3,3)/6;
vs(4) = v(3,2)/6;

% draw interpolated velocities, u*, v*

plot([x(2)-us(1)/umax*dx*vscal,x(2)+us(1)/umax*dx*vscal],[y(2)+dy/8,y(2)+dy/8],'g-')
plot([x(2)-us(2)/umax*dx*vscal,x(2)+us(2)/umax*dx*vscal],[y(3)-dy/8,y(3)-dy/8],'g-')
plot([x(3)-us(3)/umax*dx*vscal,x(3)+us(3)/umax*dx*vscal],[y(3)-dy/8,y(3)-dy/8],'g-')
plot([x(3)-us(4)/umax*dx*vscal,x(3)+us(4)/umax*dx*vscal],[y(2)+dy/8,y(2)+dy/8],'g-')

plot([x(2)+dx/8,x(2)+dx/8],[y(2)-vs(1)/vmax*dy*vscal,y(2)+vs(1)/vmax*dy*vscal],'c-')
plot([x(2)+dx/8,x(2)+dx/8],[y(3)-vs(2)/vmax*dy*vscal,y(3)+vs(2)/vmax*dy*vscal],'c-')
plot([x(3)-dx/8,x(3)-dx/8],[y(3)-vs(3)/vmax*dy*vscal,y(3)+vs(3)/vmax*dy*vscal],'c-')
plot([x(3)-dx/8,x(3)-dx/8],[y(2)-vs(4)/vmax*dy*vscal,y(2)+vs(4)/vmax*dy*vscal],'c-')

% cutcell geometric properties

A1 = 0.25*dx;
V1 = 0.5*A1^2;

V2 = dx*dy - V1;

deta = 0.5*(sqrt(V1)+sqrt(V2)); % length scale used to approximate pressure gradient on cutface

% compute divergence for each cutcell

ds(1) = ( -us(1) -vs(1) )*A1/V1;
ds(2) = ( -us(2) +vs(2) )*A1/V1;
ds(3) = (  us(3) +vs(3) )*A1/V1;
ds(4) = (  us(4) -vs(4) )*A1/V1;

ds(5) = ( (u(3,1) - u(2,1))*dy + vs(1)*A1 + vs(4)*A1 - v(2,1)*dx )/V2;
ds(6) = ( us(1)*A1 + us(2)*A1 - u(1,2)*dy + (v(1,3) - v(1,2))*dx )/V2;
ds(7) = ( (u(3,3) - u(2,3))*dy + v(2,4)*dx - vs(2)*A1 - vs(3)*A1 )/V2;
ds(8) = ( u(4,2)*dy - us(3)*A1 - us(4)*A1 + (v(3,3) - v(3,2))*dx )/V2;

% build A matrix

n_cells = 8;
A = sparse(n_cells,n_cells)














































