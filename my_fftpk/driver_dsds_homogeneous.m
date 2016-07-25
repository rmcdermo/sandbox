% McDermott
% 7-25-2016
% driver_dsds_homogeneous.m
%
% DS-DS problem with variable source term and
% homogeneous Dirichlet boundaries
% ---------------------------------------------------

% define the source term

fc = -sin(xc);

% step 1: analysis (determine fbar from inverse)

fbar = fft_dsds(fc);

% step 2: solve

ubar = solve_dsds(fbar);

% step 3: synthesis

uc = ifft_dsds(ubar) * dx^2;

plot(xc,uc,'o')

% check discrete solution of Poisson

for ii=1:nx

    % apply boundary conditions
    if ii==1
        uc_im1 = -uc(1);
    else
        uc_im1 = uc(ii-1);
    end
    if ii==nx
        uc_ip1 = -uc(nx);
    else
        uc_ip1 = uc(ii+1);
    end

    rc(ii) = ( uc_im1 - 2*uc(ii) + uc_ip1 )/dx^2 - fc(ii);
end
max(abs(rc))