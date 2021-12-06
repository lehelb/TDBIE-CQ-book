function u = evalCQ(g,Ks,N,T,dlt)
%function u = evalCQ(g,Ks,N,T,dlt)
%evaluate u_n = \sum_{j = 0}^N \omega_{n-j}(K) g(t_j) = g(t_n), 
%to solve a discrete convolution equation call with K^{-1}
if (nargin < 5)
    dlt = @(z) 2*(1-z)./(1+z); %trapezoid rule
end
dt = T/N; lam = 10^(-14/(2*N+1));
ts = (0:N)*dt; lams = lam.^(0:N);
zs = lam*exp(-2*pi*1i*(0:N)/(N+1));
u = (1./lams).*ifft(Ks(dlt(zs)/dt).*fft(lams.*g(ts)));
