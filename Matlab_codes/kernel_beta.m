function [us,us_ex,M] = kernel_beta(Q,beta,T,N)
%function [us,us_ex,M,ts] = kernel_beta(Q,beta,T,N)

n0 = 5; Ks = @(s) exp(-s.^(beta/2));
phi_lim = max(0,pi*(1-1/beta));
phi = (phi_lim+pi/2)/2;
mu = -beta/2; tol = 1e-8;
g = @(t) exp(-400*(t-0.5).^2);

[us,M] = fast_and_obl(g,N,T,Ks,mu,phi,n0,Q,tol);
us_ex = real(evalCQ(g,Ks,N,T,@(zt) 1-zt+.5*(1-zt).^2));