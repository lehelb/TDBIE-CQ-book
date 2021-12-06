function u = evalRKCQ(g,Ks,N,T,RK)
%function u = evalRKCQ(g,Ks,N,T,RK)
%RK = 0 for 2 stage Radau IIA
%RK = 1 for 3 stage Radau IIA
%RK = 2 for 3 stage Lobatto IIIC
%RK = 3 for 2 stage Gauss method

dt = T/N; lam = 10^(-14/(2*N+1));
ts = (0:N)*dt;
[A,b,c] = RKdata(RK);
m = length(c);
gs = zeros(m,N+1);
for j = 1:m
  gs(j,:) = (lam.^(0:N)).*g(ts+dt*c(j));  
end
rhss = fft(gs.').'; u_hat = zeros(1,N+1);
em0 = (A'\b)'; em1 = 1-em0*ones(m,1);
gamma = @(zt) (zt./(1-zt*em1))*em0;  
z = exp(-1i*2*pi/(N+1));zs = z.^(0:N);
for L = 1:N+1
    zt = lam*zs(L);
    [V,D] = eig((zt./(1-zt))*ones(m,1)*b'+A);
    D = diag(1./diag(D))/dt;    
    ev = Ks(diag(D));            
    u_hat_tmp = V*diag(ev)*(V\rhss(:,L));
    u_hat(L) = gamma(zt)*u_hat_tmp;    
end
u = ifft(u_hat);u = real((lam.^(-(0:N))).*u(1:N+1));
