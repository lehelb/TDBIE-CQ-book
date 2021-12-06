%%
N = 10; T = 2; RK = 3; 
g = @(t) (sin(t)).^2; 
%u_ex = @(t) 2*cos(t).*(sin(t));
Ks = @(s) s.^(-1.231);
u = evalRKCQ(g,Ks,N,T,RK);
u1 = evalRKCQ(g,Ks,N,T,1);
ts = (0:N)*T/N;plot(ts,u-u1);

%%
N = 40; T = 4; ts = (0:N)*T/N;
g = @(t) exp(-100*(t-0.5).^2);
Ks = @(s) (1./s)./(1-exp(-s));
u = evalRKCQ(g,Ks,N,T,3);
plot(ts,real(u),ts,g(ts)+g(ts-1)+g(ts-2)+g(ts-3))
legend('computed','exact')