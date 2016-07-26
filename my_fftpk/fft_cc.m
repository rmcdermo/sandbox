% McDermott
% 7-22-2016
% fft_cc.m
%
% forward transform (analysis) for cyclic boundary conditions

function [xbar] = fft_cc(x)

n = length(x);
if mod(n,2)==0
    m = n/2 - 1; % n is even
    n_even = true;
else
    m = (n-1)/2; % n is odd
    n_even = false;
end

t = 2*pi/n;

xbar = zeros(1,n);
xbar(1) = sum(x);

for jj=1:m
    xbar(2*jj)   = 0;
    xbar(2*jj+1) = 0;
    for ii=1:n
        xbar(2*jj)   = xbar(2*jj)   + x(ii) * cos(t*ii*jj);
        xbar(2*jj+1) = xbar(2*jj+1) + x(ii) * sin(t*ii*jj);
    end
end

if n_even
    for ii=1:n
        xbar(n) = xbar(n) + x(ii)*(-1)^ii;
    end
end

xbar = 2/n * xbar;