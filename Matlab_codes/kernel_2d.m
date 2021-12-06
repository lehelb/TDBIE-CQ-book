function [us,us_ex,M] = kernel_2d(Q,r,T,N)
%function [us,us_ex,M,ts] = kernel_damped(Q,r,T,N)


Ks = @(s) besselk(0,s*r);
Ks_t = @(s) besselk(0,s*r,1); %K_0 scaled by e^{sr}

phi = pi/4; mu = 0; tol = 1e-10;
g = @(t) exp(-400*(t-0.5).^2);

[us_t,M] = fast_and_obl(g,N,T,Ks_t,mu,phi,10,Q,tol);
%shift using CQ but compute only first gamma *r/dt weights
gamma = 3; dt = T/N;
N0 = ceil(gamma*r/dt);

ws = real(weightsCQ(N0,N0*dt,@(s) exp(-s*r),@(zt) 1-zt+.5*(1-zt).^2));
us = zeros(1,N+1);
for n = 0:N
    for j = 0:min(N0,n)       
        us(n+1) = us(n+1)+ws(j+1)*us_t(n-j+1);
    end
end

us_ex = real(evalCQ(g,Ks,N,T,@(zt) 1-zt+.5*(1-zt).^2));