% McDermott
% 7-25-2016
% ifft_nsns.m
%
% backward transform (synthesis) for NS-NS bcs

function [x] = ifft_nsns(xbar)

n = length(xbar);

t = pi/(2*n) * [0:n-1];

x = zeros(1,n);

for ii=1:n
    x(ii) = sum( xbar .* cos( (2*ii-1)*t ) );
end