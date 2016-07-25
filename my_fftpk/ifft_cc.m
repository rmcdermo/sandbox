% McDermott
% 7-22-2016
% ifft_cc.m
%
% backward transform (synthesis) for cyclic boundary conditions

function [x] = ifft_cc(xbar)

n = length(xbar);
if mod(n,2)==0
    m = n/2 - 1; % n is even
    n_even = 1;
else
    m = (n-1)/2; % n is odd
    n_even = 0;
end

t = 2*pi/n;
x = zeros(1,n);

if n_even

    hxb1 = 0.5*xbar(1);
    hxbn = 0.5*xbar(n);
    for ii=1:n
        sumii = 0;
        for jj=1:m
            sumii = sumii + xbar(2*jj)*cos(t*ii*jj) + xbar(2*jj+1)*sin(t*ii*jj);
        end
        x(ii) = hxb1 + sumii + hxbn*(-1)^ii;
    end

else

    hxb1 = 0.5*xbar(1);
    for ii=1:n
        sumii = 0;
        for jj=1:m
            sumii = sumii + xbar(2*jj)*cos(t*ii*jj) + xbar(2*jj+1)*sin(t*ii*jj);
        end
        x(ii) = hxb1 + sumii;
    end

end
