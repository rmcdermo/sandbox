% McDermott
% 7-25-2016
% fft_dsds.m
%
% forward transform (analysis) for DS-DS boundaries (Dirichlet Staggered)

function [xbar] = fft_dsds(x)

n = length(x);

t = (2*[1:n]-1)*pi/(2*n);

xbar = zeros(1,n);

for jj=1:n
    xbar(jj) = sum( x .* sin(t*jj) );
end

xbar = 2/n * xbar;