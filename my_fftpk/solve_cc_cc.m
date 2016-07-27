% McDermott
% 7-26-2016
% solve_cc_cc.m
%
% solve eigenvalue system for cyclic boundary conditions in 2D

function [xbar] = solve_cc_cc(ybar,dd)

n = size(ybar)

if mod(n(1),2)==0
    m(1) = n(1)/2 - 1; % n(1) is even
    n_even(1) = true;
else
    m(1) = (n(1)-1)/2; % n(1) is odd
    n_even(1) = false;
end

if mod(n(2),2)==0
    m(2) = n(2)/2 - 1; % n(2) is even
    n_even(2) = true;
else
    m(2) = (n(2)-1)/2; % n(2) is odd
    n_even(2) = false;
end

t = pi/n;

xbar(1) = 0;
for jj=1:m
    lam_1 = -4 * sin(jj*t)^2 / dd(1);
    if lam_1==0
        xbar(2*jj)   = 0;
        xbar(2*jj+1) = 0;
    else
        xbar(2*jj)   = ybar(2*jj) / lambda;
        xbar(2*jj+1) = ybar(2*jj+1) / lambda;
    end
end
if n_even
    lambda = -4 / dxdx;
    xbar(n) = ybar(n) / lambda;
end
