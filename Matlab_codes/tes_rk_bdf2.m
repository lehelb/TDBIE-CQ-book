%%
T = 4;
g = @(t) exp(-100*(t-0.5).^2);
ex = @(t) g(t)+g(t-1)+g(t-2)+g(t-3)+g(t-4)+g(t-5)+g(t-6)+g(t-7);
Ks = @(s) 1./(1-exp(-s));

Ns = [10 20 40 80 160 320];

err_RK= [];
err_bdf2 = [];

for N = Ns
    ts = (0:N)*T/N;
    u_rk = evalRKCQ(g,Ks,N,T,1); %3-stage Radau IIA method
    %u = evalCQ(g,Ks,N,T,@(zt) 1-zt+.5*(1-zt).^2);
    u = evalRKCQ(g,Ks,N,T,0);
    err_RK(end+1) = max(abs(u_rk.'-ex(ts)));
    err_bdf2(end+1) = max(abs(u.'-ex(ts)));
end

loglog(Ns*3,err_RK,Ns*2,err_bdf2)
legend('RK','BDF2')
