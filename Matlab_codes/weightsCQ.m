function w = weightsCQ(N,T, Ks, dlt)
%function w = weightsCQ(N,T, Ks, dlt)
%returns conv. weights \omega_j(K) for j = 0, .., N
%based on the trapezoidal rule unless the generating
%function dlt of the scheme is given

dt = T/N;lam = 10^(-14/(2*N+1));
if (nargin < 4)
    dlt = @(z) 2*(1-z)./(1+z); %trapezoidal rule
end
zs = exp(-1i*2*pi*(0:N)/(N+1));
w = (lam.^(0:-1:-N)).*ifft(Ks(dlt(lam*zs)/dt));
