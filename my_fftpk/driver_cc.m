% McDermott
% 7-25-2016
% driver_cc.m
%
% periodic problem with variable source term
% note: f must be periodic on (0,Lx) for this to work
% ---------------------------------------------------

% define the source term

fc = -sin(xc); % implies solution of Poisson (u'' = f) is u = sin(x)

% step 1: analysis (determine fbar from inverse)

fbar = fft_cc(fc);

% step 2: solve

ubar = solve_cc(fbar,dx^2);

% step 3: synthesis

uc = ifft_cc(ubar);

plot(xc,uc,'o')

% check discrete solution of Poisson

for ii=1:nx
    im1=ii-1;
    ip1=ii+1;
    % periodic bcs
    if im1<1;  im1=im1+nx; end
    if ip1>nx; ip1=ip1-nx; end

    rc(ii) = ( uc(im1) - 2*uc(ii) + uc(ip1) )/dx^2 - fc(ii);
end
max(abs(rc))