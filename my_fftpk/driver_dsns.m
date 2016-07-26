% McDermott
% 7-26-2016
% driver_dsns.m
%
% DS-NS problem inhomogeneous boundaries
% ---------------------------------------------------

% define the source term

fc = -2*ones(1,nx);
bxs = 5;
bxf = -10;

y = fc;
y(1)  = fc(1) - 2*bxs/dx^2; % DS boundary
y(nx) = fc(nx) - bxf/dx;    % NS boundary

% step 1: analysis (determine fbar from inverse)

fbar = fft_dsns(y);

% step 2: solve

ubar = solve_dsns(fbar);

% step 3: synthesis

uc = ifft_dsns(ubar) * dx^2;

plot(xc,uc,'o')

% check discrete solution of Poisson

for ii=1:nx

    % apply boundary conditions
    if ii==1
        uc_im1 = 2*bxs - uc(1);
    else
        uc_im1 = uc(ii-1);
    end
    if ii==nx
        uc_ip1 = uc(nx) + bxf*dx;
    else
        uc_ip1 = uc(ii+1);
    end

    rc(ii) = ( uc_im1 - 2*uc(ii) + uc_ip1 )/dx^2 - fc(ii);
end
max(abs(rc))



