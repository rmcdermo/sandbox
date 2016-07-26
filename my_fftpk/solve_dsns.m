% McDermott
% 7-26-2016
% solve_dsns.m
%
% solve eigenvalue system for DS-NS (or NS-DS) bcs

function [xbar] = solve_dsns(ybar)

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