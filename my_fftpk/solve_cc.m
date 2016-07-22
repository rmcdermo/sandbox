% McDermott
% 7-22-2016
% solve_cc.m
%
% solve eigenvalue system for cyclic boundary conditions

function [xbar] = solve_cc(ybar)

n = length(ybar);
if mod(n,2)==0
    m = n/2. - 1; % n is even
    n_even = 1;
else
    m = (n-1)/2.; % n is odd
    n_even = 0;
end
t = pi/n;
lambda = zeros(1,n);
for jj=1:m
    lambda(2*jj)   = sin(t*jj)^2;
    lambda(2*jj+1) = lambda(2*jj);
end
if n_even
    for ii=1:n
        lambda(n) = 1.;
    end
end
lambda = -4. * lambda;

% solve
xbar(2:n) = ybar(2:n)./lambda(2:n);
xbar(1) = 0;