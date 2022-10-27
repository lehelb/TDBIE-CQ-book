function [s,w] = fob_weights(Q,mu,phi,T, n0, dt, tol)
% function [s,w] = fob_weights(Q,phi,T, n0, dt, tol)
% nodes and weights for fast and oblivious quadrature
% base on BDF2

%truncate to L = A/dt with error within the tolerance tol
tr_bnd = @(A) ((dt./A).^mu/(n0*cos(phi)*pi)).*(sqrt(2*A*cos(phi)-1)-2).^(-n0);
A = fzero(@(A) tr_bnd(A)-tol, n0*5);%initial guess may need to be played with
L = A/dt;

%Gauss weights and nodes
[xg,wg] = gauss(Q); xg = (xg+1)/2; wg = wg/2;
B = 2; L0 = .5/T;
J = ceil(log(L/L0)/log(1+B))
s = zeros((J+1)*Q,1);
w = zeros(1,(J+1)*Q);

% quadrature in the (upper) circular contour
s(1:Q) = L0*exp(1i*xg*(pi-phi));
w(1:Q) = (dt*(pi-phi)/(2*pi))*L0*exp(1i*xg.'*(pi-phi)).*wg;
Ls = L0*(1+B).^(0:J);
% quadrature in the (upper) straight contour
for j = 1:J
    s((j*Q+1):(j+1)*Q) = ((Ls(j+1)-Ls(j))*xg+Ls(j))*exp(1i*(pi-phi));
    w((j*Q+1):(j+1)*Q) = (exp(1i*(pi-phi))*dt*(Ls(j+1)-Ls(j))/(2*pi*1i))*wg;
end