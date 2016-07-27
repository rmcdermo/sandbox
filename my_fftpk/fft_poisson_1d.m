% McDermott
% 7-22-2016
% fft_poisson_solver.m
%
% This code is a test bed for learning how to solve Poisson using FFT following
%
% U. Schumann and R. Sweet. Fast Fourier Transforms for Direct Solution of Poisson's
% Equation with Staggered Boundary Conditions.  J. Comput. Phys., 75:123-137, 1988.
%
% Our goal is to develop a direct solver for mixed inhomogeneous bcs.
%
% This routine selects from a set of drivers for simple 1D problems.

close all
clear all

Lx = 2*pi;
nx = 32;
dx = Lx/nx;
xf = [0:dx:Lx];
xc = xf(1:nx) + 0.5*dx;

test_case = 5

switch test_case
    case 1; driver_cc
    case 2; driver_dsds
    case 3; driver_nsns
    case 4; driver_dsns
    case 5; driver_nsds
end

