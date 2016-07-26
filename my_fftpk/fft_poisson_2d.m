% McDermott
% 7-26-2016
% fft_poisson_2d.m
%
% U. Schumann and R. Sweet. Fast Fourier Transforms for Direct Solution of Poisson's
% Equation with Staggered Boundary Conditions.  J. Comput. Phys., 75:123-137, 1988.
%
%  u(i-1,j) - 2*u(i,j) + u(i+1,j)     u(i,j-1) - 2*u(i,j) + u(i,j+1)
%  ------------------------------  +  ------------------------------  =  f(i,j)
%              dx^2                               dy^2

close all
clear all

Lx = 4*pi;
Ly = 2*pi;
nx = 64;
ny = 32;
dx = Lx/nx;
dy = Ly/ny;
xf = [0:dx:Lx];
yf = [0:dy:Ly];
xc = xf(1:nx) + 0.5*dx;
yc = yf(1:ny) + 0.5*dy;

% periodic problem with variable source term
% note: f must be periodic on (0,Lx) and (0,Ly) for this to work
% ---------------------------------------------------------------

% define the source term

for ii=1:nx
    for jj=1:ny
        fc(ii,jj) = -( sin(xc(ii)) + sin(yc(jj)) ); % implies solution of Poisson (u'' = f) is u(x,y) = sin(x) + sin(y)
    end
end

% step 1: analysis (determine fbar from inverse)

fbar=zeros(nx,ny);
for jj=1:ny
    fbar(:,jj) = fft_cc(fc(:,jj));
end
fbar=fbar';
for ii=1:nx
    fbar(:,ii) = fft_cc(fbar(:,ii));
end
fbar=fbar';

% step 2: solve

ubar = solve_cc_cc(fbar,[dx^2,dy^2]);

return

% step 3: synthesis

uc = ifft_cc(ubar);

plot(xc,uc,'o')

% check discrete solution of Poisson

for ii=1:nx
    im1=ii-1;
    ip1=ii+1;
    % periodic bcs
    if im1<1;  im1=im1+nx; end
    if ip1>nx; ip1=ip1-nx; end

    rc(ii) = ( uc(im1) - 2*uc(ii) + uc(ip1) )/dx^2 - fc(ii);
end
max(abs(rc))
