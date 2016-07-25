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
    x(ii) = sum( xbar .* sin( (2*ii-1)*t ) );
end