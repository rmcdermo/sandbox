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
%
% CAUTION: Transforms from Schumann/Sweet paper do not work for arbitrary Dirichlet
%          data on either side.  We use transforms from FLASH version of fftpack.f90.

close all
clear all

% define the grid

Lx = 2*pi;
nx = 32;
dx = Lx/nx;
xf = [0:dx:Lx];
xc = xf(1:nx) + 0.5*dx;

test_case = 'dirichlet inhomogeneous'

switch test_case
    case 'periodic';                 driver_cc
    case 'dirichlet homogeneous';    driver_dsds_homogeneous
    case 'dirichlet inhomogeneous';  driver_dsds_inhomogeneous
end

