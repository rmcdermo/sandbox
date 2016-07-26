% McDermott
% 7-26-2016
% driver_nsns.m
%
% NS-NS problem inhomogeneous boundaries
% ---------------------------------------------------

% define the source term

fc = zeros(1,nx);
bxs = 2;
bxf = -3;

y = fc;
y(1)  = fc(1)  + bxs/dx;
y(nx) = fc(nx) - bxf/dx;

pois_ptb = mean(y);

% step 1: analysis (determine fbar from inverse)

fbar = fft_nsns(y);

% step 2: solve

ubar = solve_nsns(fbar,dx^2);

% step 3: synthesis

uc = ifft_nsns(ubar);

plot(xc,uc,'o')

% check discrete solution of Poisson

for ii=1:nx

    % apply boundary conditions
    if ii==1
        uc_im1 = uc(1)  - bxs*dx;
    else
        uc_im1 = uc(ii-1);
    end
    if ii==nx
        uc_ip1 = uc(nx) + bxf*dx;
    else
        uc_ip1 = uc(ii+1);
    end

    rc(ii) = ( uc_im1 - 2*uc(ii) + uc_ip1 )/dx^2 - fc(ii) + pois_ptb;
end
max(abs(rc))



