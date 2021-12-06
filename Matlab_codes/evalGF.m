function u = evalGF(gs,Kzt)
%function u = evalGF(g,Kzt)
%evaluate u_n = \sum_{j = 0}^N \omega_{n-j} g_j 
% where \omega_j are given by the generating function
% Kzt(\zeta). 

N = length(gs)-1;
lam = 10^(-14/(2*N+1)); lams = lam.^(0:N);
zs = lam*exp(-2*pi*1i*(0:N)/(N+1));
u = (1./lams).*ifft(Kzt(zs).*fft(lams.*gs));
