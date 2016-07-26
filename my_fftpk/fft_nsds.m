% McDermott
% 7-26-2016
% fft_dsns.m
%
% forward transform (analysis) for NS-DS boundaries
% see Schumann and Sweet Eq. (25)

function [xbar] = fft_nsds(x)

n = length(x);

t = pi/(4*n);

xbar = zeros(1,n);

for jj=1:n
    xbar(jj) = 0;
    for ii=1:n
        xbar(jj) = xbar(jj) + x(ii) * cos((2*ii-1)*(2*jj-1)*t);
    end
end

xbar = 2/n * xbar;

