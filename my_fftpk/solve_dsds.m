% McDermott
% 7-25-2016
% solve_dsds.m
%
% solve eigenvalue system for DS-DS bcs

function [xbar] = solve_dsds(ybar)

n = length(ybar);
t = pi/(2*n);
for jj=1:n
    lambda = -4 * sin(jj*t)^2; % compute eigenvalue, SS Eq. (32)
    if (lambda==0)
        xbar(jj) = 0;
    else
        xbar(jj) = ybar(jj) / lambda; % solve
    end
end