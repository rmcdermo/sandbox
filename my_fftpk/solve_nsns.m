% McDermott
% 7-25-2016
% solve_nsns.m
%
% solve eigenvalue system for NS-NS bcs

function [xbar] = solve_nsns(ybar)

n = length(ybar);
t = pi/(2*n);
for jj=1:n
    lambda = -4 * sin((jj-1)*t)^2; % compute eigenvalue, SS Eq. (33)
    if (lambda==0)
        xbar(jj) = 0;
    else
        xbar(jj) = ybar(jj) / lambda; % solve
    end
end