% McDermott
% 7-22-2016
% fft_poisson_solver.m
%
% This code is a test bed for learning how to solve Poisson and Laplace using FFT.
%
% U. Schumann and R. Sweet. Fast Fourier Transforms for Direct Solution of Poisson's
% Equation with Staggered Boundary Conditions.  J. Comput. Phys., 75:123-137, 1988.
%
% Consider the Poisson equation
%
%   u_{i-1} - 2u_i + u_{i+1} = f_i
%
% We will look at the case where y_i is non-zero (Poisson) and y_i = 0 (Laplace).
% We will also look at *inhomogeneous* Dirichlet and Neumann boundaries on a
% staggered grid, DS and NS boundaries in the Schumann-Sweet paper.

close all
clear all

% define the grid

Lx = 2*pi;
nx = 32;
dx = Lx/nx;
xf = [0:dx:Lx];
xc = xf(1:nx) + 0.5*dx;


% % periodic problem with variable source term
% % note: f must be periodic on (0,Lx) for this to work
% % ---------------------------------------------------

% % define the source term

% fc = -sin(xc); % implies solution of Poisson (u'' = f) is u = sin(x)

% % step 1: analysis (determine fbar from inverse)

% fbar = fft_cc(fc);

% % step 2: solve

% ubar = solve_cc(fbar);

% % step 3: synthesis

% uc = ifft_cc(ubar) * dx^2;

% plot(xc,uc,'o')

% % check discrete solution of Poisson

% for ii=1:nx
%     im1=ii-1;
%     ip1=ii+1;
%     % periodic bcs
%     if im1<1;  im1=im1+nx; end
%     if ip1>nx; ip1=ip1-nx; end

%     rc(ii) = ( uc(im1) - 2*uc(ii) + uc(ip1) )/dx^2 - fc(ii);
% end
% max(abs(rc))


% DS-DS problem with variable source term
% ---------------------------------------------------

% define the source term

fc = -sin(xc);
bxs = 0;
bxf = 0;

% step 1: analysis (determine fbar from inverse)

fbar = fft_dsds(fc);

% step 2: solve

ubar = solve_dsds(fbar);

% step 3: synthesis

uc = ifft_dsds(ubar) * dx^2;

plot(xc,uc,'o')

% check discrete solution of Poisson

for ii=1:nx

    % apply boundary conditions
    if ii==1
        uc_im1 = 2*bxs - uc(1);
    else
        uc_im1 = uc(ii-1);
    end
    if ii==nx
        uc_ip1 = 2*bxf - uc(nx);
    else
        uc_ip1 = uc(ii+1);
    end

    rc(ii) = ( uc_im1 - 2*uc(ii) + uc_ip1 )/dx^2 - fc(ii);
end
max(abs(rc))












