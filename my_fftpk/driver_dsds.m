% McDermott
% 7-25-2016
% driver_dsds.m
%
% DS-DS problem with variable source term and
% non-zero Dirichlet boundaries
% ---------------------------------------------------

% define the source term

fc = -2*ones(1,nx);
%fc = -sin(xc);
%fc = zeros(1,nx);
bxs = -3;
bxf = 1;

y = fc;
y(1)  = fc(1)  - 2*bxs/dx^2;
y(nx) = fc(nx) - 2*bxf/dx^2;

% step 1: analysis (determine fbar from inverse)

fbar = fft_dsds(y);

% step 2: solve

ubar = solve_dsds(fbar);

% step 3: synthesis

uc = ifft_dsds(ubar) * dx^2;

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
        uc_ip1 = 2*bxf - uc(nx);
    else
        uc_ip1 = uc(ii+1);
    end

    rc(ii) = ( uc_im1 - 2*uc(ii) + uc_ip1 )/dx^2 - fc(ii);
end
max(abs(rc))



