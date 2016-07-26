% McDermott
% 7-22-2016
% solve_cc.m
%
% solve eigenvalue system for cyclic boundary conditions

function [xbar] = solve_cc(ybar,dxdx)

n = length(ybar);
if mod(n,2)==0
    m = n/2 - 1; % n is even
    n_even = 1;
else
    m = (n-1)/2; % n is odd
    n_even = 0;
end

t = pi/n;

xbar(1)   = 0;

for jj=1:m
    lambda = -4 * sin(jj*t)^2 / dxdx;
    if lambda==0
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
