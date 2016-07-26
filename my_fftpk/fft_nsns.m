% McDermott
% 7-25-2016
% fft_nsns.m
%
% forward transform (analysis) for NS-NS boundaries (Neumann Staggered)

function [xbar] = fft_nsns(x)

n = length(x);

t = pi/(2*n);

xbar = zeros(1,n);

for jj=1:n
    xbar(jj) = 0;
    for ii=1:n
        xbar(jj) = xbar(jj) + x(ii) * cos( (2*ii-1)*(jj-1)*t );
    end
end

xbar(1)   = 1/n * xbar(1);
xbar(2:n) = 2/n * xbar(2:n);