% McDermott
% 7-25-2016
% ifft_dsds.m
%
% backward transform (synthesis) for DS-DS bcs
% taken from Schumann and Sweet Eq. (11)
% note: error in last term in SS

function [x] = ifft_dsds(xbar)

n = length(xbar);

t = pi/(2*n);

x = zeros(1,n);

for ii=1:n
    x(ii) = 0;
    for jj=1:n-1
        x(ii) = x(ii) + xbar(jj) * sin( (2*ii-1)*t*jj );
    end
    x(ii) = x(ii) + (-1)^(ii-1) * xbar(n) / 2;
end