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
xf = 0:dx:Lx;
xc = xf(1:nx) + 0.5*dx;
x  = 2*pi*xc/Lx;

% define the source term
f = -sin(x); % implies solution of Poisson (u'' = f) is u = sin(x)
y = dx^2 * f;

% periodic problem with variable source term
% ------------------------------------------

% step 1: analysis (determine fbar from inverse)

ybar = fft_cc(y);

% fvar=fft(f);

% fbar2(1)=real(fvar(1));

% for i=1:nx/2
%     fbar2(2*i)  =real(fvar(i+1));
%     fbar2(2*i+1)=imag(fvar(i+1));
% end

% return

% step 2: solve

ubar = solve_cc(ybar);

% step 3: synthesis

u = ifft_cc(ubar);

plot(x,u,'o')