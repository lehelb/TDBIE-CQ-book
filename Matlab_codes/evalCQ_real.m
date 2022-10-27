function u = evalCQ_real(g,Ks,N,T,dlt)
%function u = evalCQ_real(g,Ks,N,T,dlt)
%evaluate u_n = \sum_{j = 0}^N \omega_{n-j}(K) g(t_j) = g(t_n)
%in the case of real u and g.
%To solve a discrete convolution equation call with K^{-1}
if (nargin < 5)
    dlt = @(z) 2*(1-z)./(1+z); %trapezoidal rule
end
dt = T/N; lam = 10^(-14/(2*N+1));
ts = (0:N)*dt; lams = lam.^(0:N);
Nhalf = floor(N/2);
zs = lam*exp(-2*pi*1i*(0:Nhalf)/(N+1));
Ks_eval = zeros(1,N+1);
Ks_eval(1:Nhalf+1) = Ks(dlt(zs)/dt);% compute Ks for half of the freq.
Ks_eval(N+1:-1:Nhalf+2) = conj(Ks_eval(2:N-Nhalf+1)); %the other half are conjugates
u = (1./lams).*real(ifft(Ks_eval.*fft(lams.*g(ts))));
