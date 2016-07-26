% McDermott
% 7-26-2016
% ifft_dsns.m
%
% backward transform (synthesis) for DS-NS bcs
% see Schumann and Sweet (JCP, 1988) Eq. (13)

function [x] = ifft_dsns(xbar)

n = length(xbar);

t = pi/(4*n);

x = zeros(1,n);

for ii=1:n
    x(ii) = 0;
    for jj=1:n
        x(ii) = x(ii) + xbar(jj) * sin( (2*ii-1)*(2*jj-1)*t );
    end
end