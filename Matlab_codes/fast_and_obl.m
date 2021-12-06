function [us,M] = fast_and_obl(g,N,T,Ks,mu,phi,n0,Q,tol)
%function [us,M] = fast_and_obl(g,N,T,Ks,mu,phi,n0,Q,tol)
% BDF2 fast and oblivious quadrature computation of
% K(\partial_t) g

dt = T/N;
dlt = @(zt) 1-zt+.5*(1-zt).^2; %BDF2 CQ
dw = (1/dt)*[3/2 -2 .5]; %coeff of BDF2
ej = @(z,j) ((2-sqrt(1+2*z)).^(-j-1)-(2+sqrt(1+2*z)).^(-j-1))./sqrt(1+2*z);
ws(1:n0) = real(weightsCQ(n0-1,(n0-1)*dt, Ks,dlt));
[S,W] = fob_weights(Q,mu,phi,T, dt, @(z) ej(z,n0-1),tol);

M = length(S);
ys = zeros(M,2);us = zeros(1,N+1);

for n = 0:n0-1
    for j = 0:n
        us(n+1) = us(n+1)+ws(j+1)*g(dt*(n-j));        	
    end    
end

for n = 0:n0
    tmp = ys(:,2);
    ys(:,2) = (1./(dw(1)-S)).*(-dw(2)*ys(:,2)-dw(3)*ys(:,1)+g(dt*n));   
    ys(:,1) = tmp;
end

for n = n0:N    
    corr = zeros(M,1);
    for j = 0:(n0-1)
        us(n+1) = us(n+1)+ws(j+1)*g(dt*(n-j));        
        corr = corr+ej(S*dt,j)*g(dt*(n-j));
    end    
    us(n+1) = us(n+1) + 2*real(W*(Ks(S).*(ys(:,2)/dt-corr)));    

    tmp = ys(:,2);    
    ys(:,2) = (1./(dw(1)-S)).*(-dw(2)*ys(:,2)-dw(3)*ys(:,1)+g(dt*(n+1)));   
    ys(:,1) = tmp;
end
