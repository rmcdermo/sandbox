% McDermott
% 7-26-2016
% solve_nsds.m
%
% solve eigenvalue system for DS-NS (or NS-DS) bcs
% Schumann and Sweet (JCP, 1988) Eq. (34)

function [xbar] = solve_nsds(ybar)

n = length(ybar);
t = pi/(4*n);
for jj=1:n
    lambda = -4 * sin((2*jj-1)*t)^2;
    if (lambda==0)
        xbar(jj) = 0;
    else
        xbar(jj) = ybar(jj) / lambda;
    end
end