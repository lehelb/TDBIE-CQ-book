function x = solveLT(N1,N2,g,Ks,dt,rhs,dlt)
%function x = solveLT(N1,N2,g,Ks,dt,rhs)
%solve \sum_{j = 0}^N \omega_{n-j} x_j = g(t_{N_1}+n)+rhs_n, n =
%0,\dots,N where N = N_2-N_1

if (nargin < 7)
    dlt = @(z) 2*(1-z)./(1+z); %trapezoidal rule
end
N = N2-N1;

if (nargin < 6)
    rhs = zeros(N+1,1);
end

if (N == 0)    
    x = Ks(dlt(0)/dt)\(g(N1*dt)+rhs(1));
else
    Nhalf = floor((N1+N2)/2);            
    x1 = solveLT(N1,Nhalf,g,Ks,dt,rhs(1:(Nhalf-N1+1)),dlt);
    
    lam = 10^(-14/(2*N+1));zs = dlt(lam*exp(-(0:N)'*2*pi*1i/(N+1)))/dt;
    x1_tmp = zeros(N+1,1); x1_tmp(1:length(x1)) = x1;
    rhs_n = (lam.^(0:-1:-N)').*ifft(Ks(zs).*fft(x1_tmp.*(lam.^(0:N)')));
    
    rhs_n = rhs(Nhalf-N1+2:end)-rhs_n(Nhalf-N1+2:end);    
    x2 = solveLT(Nhalf+1,N2,g,Ks,dt,rhs_n,dlt);
    x = [x1;x2];
end
