% McDermott
% 7-25-2016
% ifft_dsds.m
%
% backward transform (synthesis) for DS-DS bcs

function [x] = ifft_dsds(xbar)

n = length(xbar);

t = pi/(2*n) * [1:n];

x = zeros(1,n);

for ii=1:n
    x(ii) = (-1)^(ii-1) * xbar(n) + 2 * sum( xbar(1:n-1) .* sin( (2*ii-1)*t(1:n-1) ) );
end