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

s = 0.25*dx;
v1 = 0.5*s^2;
v2 = dx*dy - v1;

de = 0.5*(sqrt(v1)+sqrt(v2)); % length scale used to approximate pressure gradient on cutface

% compute divergence for each cutcell

b(1) = -( -us(1) -vs(1) )*s;
b(2) = -( -us(2) +vs(2) )*s;
b(3) = -(  us(3) +vs(3) )*s;
b(4) = -(  us(4) -vs(4) )*s;

b(5) = -( (u(3,1) - u(2,1))*dy + vs(1)*s + vs(4)*s - v(2,1)*dx );
b(6) = -( us(1)*s + us(2)*s - u(1,2)*dy + (v(1,3) - v(1,2))*dx );
b(7) = -( (u(3,3) - u(2,3))*dy + v(2,4)*dx - vs(2)*s - vs(3)*s );
b(8) = -( u(4,2)*dy - us(3)*s - us(4)*s + (v(3,3) - v(3,2))*dx );

% build A matrix

n_cells = 8;
A = sparse(n_cells,n_cells);

% diagonal coefficients
A(1,1) = 2*s/de;
A(2,2) = 2*s/de;
A(3,3) = 2*s/de;
A(4,4) = 2*s/de;
A(5,5) = 2*s/de;
A(6,6) = 2*s/de;
A(7,7) = 2*s/de;
A(8,8) = 2*s/de;

% connectivity
A(1,5) = -s/de;
A(1,6) = -s/de;

A(2,6) = -s/de;
A(2,7) = -s/de;

A(3,7) = -s/de;
A(3,8) = -s/de;

A(4,5) = -s/de;
A(4,8) = -s/de;

A(5,4) = -s/de;
A(5,1) = -s/de;

A(6,1) = -s/de;
A(6,2) = -s/de;

A(7,2) = -s/de;
A(7,3) = -s/de;

A(8,3) = -s/de;
A(8,4) = -s/de;

% solve Poisson equation for "condensed" system
Ac = A(1:n_cells-1,1:n_cells-1);
bc = b(1:n_cells-1);
pc = Ac\bc';
p = [pc;0];

figure(2)
spy(full(A))

% project velocities

un(1) = us(1) - (p(1)-p(6))/de;
un(2) = us(2) - (p(2)-p(6))/de;
un(3) = us(3) + (p(3)-p(8))/de;
un(4) = us(4) + (p(4)-p(8))/de;

vn(1) = vs(1) - (p(1)-p(5))/de;
vn(2) = vs(2) + (p(2)-p(7))/de;
vn(3) = vs(3) + (p(3)-p(7))/de;
vn(4) = vs(4) - (p(4)-p(5))/de;

% draw new velocities, un, vn

figure(3)

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

plot([x(2)-un(1)/umax*dx*vscal,x(2)+un(1)/umax*dx*vscal],[y(2)+dy/8,y(2)+dy/8],'g-','LineWidth',5)
plot([x(2)-un(2)/umax*dx*vscal,x(2)+un(2)/umax*dx*vscal],[y(3)-dy/8,y(3)-dy/8],'g-','LineWidth',5)
plot([x(3)-un(3)/umax*dx*vscal,x(3)+un(3)/umax*dx*vscal],[y(3)-dy/8,y(3)-dy/8],'g-','LineWidth',5)
plot([x(3)-un(4)/umax*dx*vscal,x(3)+un(4)/umax*dx*vscal],[y(2)+dy/8,y(2)+dy/8],'g-','LineWidth',5)

plot([x(2)+dx/8,x(2)+dx/8],[y(2)-vn(1)/vmax*dy*vscal,y(2)+vn(1)/vmax*dy*vscal],'c-','LineWidth',5)
plot([x(2)+dx/8,x(2)+dx/8],[y(3)-vn(2)/vmax*dy*vscal,y(3)+vn(2)/vmax*dy*vscal],'c-','LineWidth',5)
plot([x(3)-dx/8,x(3)-dx/8],[y(3)-vn(3)/vmax*dy*vscal,y(3)+vn(3)/vmax*dy*vscal],'c-','LineWidth',5)
plot([x(3)-dx/8,x(3)-dx/8],[y(2)-vn(4)/vmax*dy*vscal,y(2)+vn(4)/vmax*dy*vscal],'c-','LineWidth',5)

% compute divergence for each cutcell

bn(1) = ( -un(1) -vn(1) )*s;
bn(2) = ( -un(2) +vn(2) )*s;
bn(3) = (  un(3) +vn(3) )*s;
bn(4) = (  un(4) -vn(4) )*s;

bn(5) = ( (u(3,1) - u(2,1))*dy + vn(1)*s + vn(4)*s - v(2,1)*dx );
bn(6) = ( un(1)*s + un(2)*s - u(1,2)*dy + (v(1,3) - v(1,2))*dx );
bn(7) = ( (u(3,3) - u(2,3))*dy + v(2,4)*dx - vn(2)*s - vn(3)*s );
bn(8) = ( u(4,2)*dy - un(3)*s - un(4)*s + (v(3,3) - v(3,2))*dx );

bn











































