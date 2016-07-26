% McDermott
% 7-26-2016
% driver_nsds.m
%
% NS-DS problem inhomogeneous boundaries
% ---------------------------------------------------

% define the source term

fc = -2*ones(1,nx);
bxs = 10;
bxf = 5;

y = fc;
y(1)  = fc(1)  + bxs/dx;     % NS boundary
y(nx) = fc(nx) - 2*bxf/dx^2; % DS boundary

% step 1: analysis (determine fbar from inverse)

fbar = fft_nsds(y);

% step 2: solve

ubar = solve_nsds(fbar);

% step 3: synthesis

uc = ifft_nsds(ubar) * dx^2;

plot(xc,uc,'o')

% check discrete solution of Poisson

for ii=1:nx

    % apply boundary conditions
    if ii==1
        uc_im1 = uc(1) - bxs*dx;
    else
        uc_im1 = uc(ii-1);
    end
    if ii==nx
        uc_ip1 = 2*bxf - uc(nx);
    else
        uc_ip1 = uc(ii+1);
    end

    rc(ii) = ( uc_im1 - 2*uc(ii) + uc_ip1 )/dx^2 - fc(ii);
end
max(abs(rc))



