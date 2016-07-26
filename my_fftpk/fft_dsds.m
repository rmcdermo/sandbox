% McDermott
% 7-25-2016
% fft_dsds.m
%
% forward transform (analysis) for DS-DS boundaries (Dirichlet Staggered)
% see Schumann and Sweet Eq. (22)

function [xbar] = fft_dsds(x)

n = length(x);

t = pi/(2*n);

xbar = zeros(1,n);

for jj=1:n
    xbar(jj) = 0;
    for ii=1:n
        xbar(jj) = xbar(jj) + x(ii) * sin(t*(2*ii-1)*jj);
    end
end

xbar = 2/n * xbar;