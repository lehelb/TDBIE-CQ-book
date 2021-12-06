function [s,w] = fob_weights(Q,mu,phi,T, dt, en0, tol)
% function [s,w] = fob_weights(Q,phi,T, en0, tol)
% nodes and weights for fast and oblivious quadrature

A = 1;
while (A < 1000 && abs(en0(A*exp(1i*(pi-phi)))*(A/dt)^(mu)) > tol)
    A = A+.1;
end
L = A/dt;

[xg,wg] = gauss(Q); xg = (xg+1)/2; wg = wg/2;
B = 2; L0 = .5/T;
J = ceil(log(L/L0)/log(1+B));
s = zeros((J+1)*Q,1);
w = zeros(1,(J+1)*Q);

s(1:Q) = L0*exp(1i*xg*(pi-phi));
w(1:Q) = (dt*(pi-phi)/(2*pi))*L0*exp(1i*xg.'*(pi-phi)).*wg;
Ls = L0*(1+B).^(0:J);
for j = 1:J
    s((j*Q+1):(j+1)*Q) = ((Ls(j+1)-Ls(j))*xg+Ls(j))*exp(1i*(pi-phi));
    w((j*Q+1):(j+1)*Q) = (exp(1i*(pi-phi))*dt*(Ls(j+1)-Ls(j))/(2*pi*1i))*wg;
end