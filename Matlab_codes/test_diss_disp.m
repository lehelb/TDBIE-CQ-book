%%
set_param_cq;

Ks = @(s) 1-exp(-s);
T = 4;

js = [1 2 4 8]/2;
%js = [1/2 1];
%% BDF2
Ns = 5;

for j = js
    %10*j = 1/(sqrt(2)*sigma) --> sigma = 1/(10*j*sqrt(2))
    sigma = 1/(10*j*sqrt(2));
    g = @(t) exp(-((t-0.75).^2)/(2*sigma^2));
    ex = @(t) (g(t)+g(t-1)+g(t-2)+g(t-3));
    N = ceil(Ns(end)*2^1.3); err = inf;  
    while (N < 80200) && (err > 1e-2) 
        u = evalCQ(g,@(s) 1./Ks(s),N,T,@(zt) 1-zt+.5*(1-zt).^2);    
        ts = (0:N)*T/N;    
        err = max(abs(real(u)-ex(ts)));    
        N = N+2*j;
    end
    Ns(end+1) = N;
    err
end
Ns = Ns(2:end)
%% trapezoidal
Ns_tr = 5;

for j = js
    sigma = 1/(10*j*sqrt(2));
    g = @(t) exp(-((t-0.75).^2)/(2*sigma^2));
    ex = @(t) (g(t)+g(t-1)+g(t-2)+g(t-3));
    N = ceil(Ns_tr(end)*2^(1.3)); err = inf;  
    while (N < 72000) && (err > 1e-2) 
        u = evalCQ(g,@(s) 1./Ks(s),N,T);    
        ts = (0:N)*T/N;    
        err = max(abs(real(u)-ex(ts)));    
        N = N+2*j;
    end
    Ns_tr(end+1) = N;
    err
end
Ns_tr = Ns_tr(2:end)
%% truncated trapezoidal
cs =[0.894431935027534   0.683269582897605   0.610333312441433];
dlt = @(zt) (1-zt)+.5*(1-zt).^2+(1/4)*cs(1)*(1-zt).^3+...
    (1/8)*cs(2)*(1-zt).^4+(1/16)*cs(3)*(1-zt).^5;    

Ns_ttr = 5;

for j = js
    sigma = 1/(10*j*sqrt(2));
    g = @(t) exp(-((t-0.75).^2)/(2*sigma^2));
    ex = @(t) (g(t)+g(t-1)+g(t-2)+g(t-3));
    N = ceil(Ns_ttr(end)*2^1.33); err = inf;  
    while (N < 52000) && (err > 1e-2) 
        u = evalCQ(g,@(s) 1./Ks(s),N,T,dlt);    
        ts = (0:N)*T/N;    
        err = max(abs(real(u)-ex(ts)));    
        N = N+2*j;
    end
    Ns_ttr(end+1) = N;
    err
end
Ns_ttr = Ns_ttr(2:end)
%% 2-stage Radau IIA
NsRK_2 = 5;

for j = js
    sigma = 1/(10*j*sqrt(2));
    g = @(t) exp(-((t-0.75).^2)/(2*sigma^2));
    ex = @(t) (g(t)+g(t-1)+g(t-2)+g(t-3));
    N = NsRK_2(end)*2; err = inf;  
    while (N < 6000) && (err > 1e-2) 
        u = evalRKCQ(g,@(s) 1./Ks(s),N,T,0);    
        ts = (0:N)*T/N;    
        err = max(abs(real(u)-ex(ts)));    
        N = N+2*j;
    end
    NsRK_2(end+1) = N;
    err
end
NsRK_2 = NsRK_2(2:end)

%% 3-stage Radau IIA
NsRK = 5;

for j = js
    sigma = 1/(10*j*sqrt(2));
    g = @(t) exp(-((t-0.75).^2)/(2*sigma^2));
    ex = @(t) (g(t)+g(t-1)+g(t-2)+g(t-3));
    N = NsRK(end)*2; err = inf;  
    while (N < 3000) && (err > 1e-2) 
        u = evalRKCQ(g,@(s) 1./Ks(s),N,T,1);    
        ts = (0:N)*T/N;    
        err = max(abs(real(u)-ex(ts)));    
        N = N+2*j;
    end
    NsRK(end+1) = N;
    err
end
NsRK = NsRK(2:end)

%%
sigmas = 1./((10*js)*sqrt(2));
deltas = 1./sigmas;
loglog(deltas,Ns); hold on
loglog(deltas,Ns_tr,'--');
loglog(deltas,Ns_ttr,'.-');
loglog(deltas,NsRK_2,'-.');
loglog(deltas,NsRK,'-x');

hold off; 
h = legend('BDF2','TR','TTR','RK(2)','RK(3)');
set(h,'Interpreter','Latex')
shg
h = xlabel('$\sigma^{-1}$'); set(h,'Interpreter','Latex');  h = ylabel('$N$');set(h,'Interpreter','Latex')

%% modified CQ with r
Nsmod = 5;
js = [1/2 1 2 4];
r = 1/sqrt(2);T = 4;
dlt_bdf2 = @(zt) 1-zt+.5*(1-zt).^2;
errs1 = []; Ns1 = [];
for j = js
    sigma = 1/(10*j*sqrt(2));
    g = @(t) exp(-((t-0.75).^2)/(2*sigma^2));
    %g = @(t) (1/sqrt(2*pi*sigma^2))*exp(-((t-0.75).^2)/(2*sigma^2));
    ex = @(t) (g(t)+g(t-r)+g(t-2*r)+g(t-3*r)+g(t-4*r)+g(t-5*r));
    N = ceil(Nsmod(end)); err = inf;  
    while (N < 6000) && (err > 1e-2) 
        ts = (0:N)*T/N; dt = T/N;
        n = floor(r/dt);
        gs = g(ts);
        u = evalGF(gs,@(zt) 1./(1-(zt.^(n)).*exp(-(r-n*dt)*dlt_bdf2(zt)/dt)));    
        %u = evalCQ(g,@(s) 1./(1-exp(-s*r)),N,T,dlt_bdf2);        
        ts = (0:N)*T/N;    
        err = max(abs(real(u)-ex(ts)));
        if (j == 1)
            errs1(end+1) = err; Ns1(end+1) = N;
        end
        N = N+1;
    end
    Nsmod(end+1) = N;
    err
end
Nsmod = Nsmod(2:end)
%% RK with r
NsRk = 5;
js = [1/2 1 2 4];
r = 1/sqrt(2);T = 4;
for j = js
    sigma = 1/(10*j*sqrt(2));
    %g = @(t) (1/sqrt(2*pi*sigma^2))*exp(-((t-0.75).^2)/(2*sigma^2));
    g = @(t) exp(-((t-0.75).^2)/(2*sigma^2));
    ex = @(t) (g(t)+g(t-r)+g(t-2*r)+g(t-3*r)+g(t-4*r)+g(t-5*r));
    N = NsRk(end); err = inf;  
    while (N < 6000) && (err > 1e-2)        
        u = evalRKCQ(g,@(s) 1./(1-exp(-s*r)),N,T,1);
        ts = (0:N)*T/N;    
        err = max(abs(real(u)-ex(ts)));
        %err = (T/N)*sum(abs(u-ex(ts)));
        N = N+2;
    end
    NsRk(end+1) = N;
    err
end
NsRk = NsRk(2:end)
%%
sigmas = 1./((10*js)*sqrt(2));
deltas = 1./sigmas;
loglog(deltas,3*NsRk); hold on
loglog(deltas,Nsmod,'--');


hold off; 
h = legend('RK(3)','modified CQ');
set(h,'Interpreter','Latex')
shg
h = xlabel('$\sigma^{-1}$'); set(h,'Interpreter','Latex');  h = ylabel('$N$');set(h,'Interpreter','Latex')
%%
%Ns mod conv order
N = Nsmod(1);
j = js(end);
sigma = 1/(10*j*sqrt(2));
g = @(t) exp(-((t-0.75).^2)/(2*sigma^2));
ex = @(t) (g(t)+g(t-r)+g(t-2*r)+g(t-3*r)+g(t-4*r)+g(t-5*r));
errs1 = [];
Ns1 = [];

for l = 1:6
    N = 2*N;
    ts = (0:N)*T/N; dt = T/N;
    n = floor(r/dt);
    gs = g(ts);
    u = evalGF(gs,@(zt) 1./(1-(zt.^(n)).*exp(-(r-n*dt)*dlt_bdf2(zt)/dt)));
    errs1(end+1) = max(abs(u-ex(ts)));
    Ns1(end+1) = N;
end
loglog(Ns1,errs1);shg; hold on; loglog(Ns1,100000*Ns1.^(-3)); hold off
