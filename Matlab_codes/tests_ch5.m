set(0,'defaultfigureposition',[380 320 540 200],...
'defaultaxeslinewidth',0.9,'defaultaxesfontsize',8,...
'defaultlinelinewidth',1.1,'defaultpatchlinewidth',1.1,...
'defaultlinemarkersize',15); format compact
%%
% As an example of use of |evalRKCQ| we choose a hyperbolic symbol
%
% $$ K(s) = 1-e^{-s}$$
%
% and solve the equation
%
% $$K(\partial_t) u = g.$$
%
% As $K^{-1}(s) = \frac1{1-e^{-s}} = \sum_{j = 0}^\infty
% e^{-js}$ the exact solution is given by
%
% $$u(t) = K^{-1}(\partial_t) g = \sum_{j = 0}^\infty g(t-j).$$
%
% We will apply |evalRKCQ| with $K^{-1}(s)$ to compute the 3-stage Radau IIA 
% approximation.
N = 100; T = 4; ts = (0:N)*T/N;
g = @(t) exp(-100*(t-0.5).^2);
Ks = @(s) 1./(1-exp(-s));
u = evalRKCQ(g,Ks,N,T,1); %3-stage Radau IIA method
plot(ts,real(u),ts,g(ts)+g(ts-1)+g(ts-2)+g(ts-3),'--')
legend('computed','exact','Location','west')
%%
% We see that the approximation is much better and with a much larger $\tstep$
% then the same example with the linear multistep based CQ even when
% counting the number of stages. The damping, though present is here not
% visible.
set_param_cq
T = 4;
ex_fn = @(ts) g(ts)+g(ts-1)+g(ts-2)+g(ts-3)+g(ts-4)+g(ts-5)+g(ts-6)+g(ts-7);
Ns = 12*[2 4 8 16 32 64];
err = zeros(4,length(Ns));
for j = 1:length(Ns)
    N = Ns(j); ts = (0:N)*T/N;ts1 = (0:(N/3))*T/(N/3);ts2 = (0:(N/2))*T/(N/2);
    u1 = evalRKCQ(g,Ks,N/3,T,1); %3-stage Radau IIA method
    u2 = evalRKCQ(g,Ks,N/2,T,0); %2-stage Radau IIA method
    ubdf = evalCQ(g,Ks,N,T,@(zt) 1-zt+.5*(1-zt).^2);   
    utr = evalCQ(g,Ks,N,T);    
    ex2 = ex_fn(ts2); ex1 = ex_fn(ts1); ex = ex_fn(ts);    
    err(1,j) = max(abs(u1-ex1)); err(2,j) = max(abs(u2-ex2)); 
    err(3,j) = max(abs(ubdf-ex)); err(4,j) = max(abs(utr-ex)); 
end
lini = ["-", "--", ".-", "-."];
for j = 1:4
    loglog(Ns,err(j,:),lini(j)); hold on
end
hold off; shg
xlabel('$mN$','interpreter','Latex'); ylabel('max error','interpreter','Latex')
legend('3-stage Radau IIA','2-stage Radau IIA','Trapezoidal','BDF2')
saveas(gcf,'html/conv_RK_MS','epsc')
%%
% Let us now illustrate Theorem~\ref{}. Setting $K_1(s) =
% \frac{s}{1-e^{-s}}$, we have 
% $$u_1(t) = K_1(\partial_t) g = \sum_{j = 0}^\infty g'(t-j).$$
% Making use of the 3-stage Radau IIA method, which has stage order $q = 3$ and
% full order $p = 5$ we expect convergence order $O(\Delta t^{3})$. 
% Let us also consider $K_2(s) = 1+exp(-s)$, with the exact result
% given by $u_2(t) = K(\partial_t)g(t) = g(t)+g(t-1)$. As explained in the
% previous section, we expect to see optimal convergence order.
gd = @(t) -200*(t-0.5).*g(t);
 Ks1 = @(s) s./(1-exp(-s));Ks2 = @(s) 1+exp(-s);
Ns = 2.^(1:10); T = 2;
err1 = zeros(10,1);err2 = err1; 
for j = 1:10
    N = Ns(j);ts = (0:N)*T/N;    
    u1 = evalRKCQ(g,Ks1,N,T,1);
    u2 = evalRKCQ(g,Ks2,N,T,1);
    err1(j) = norm(u1-(gd(ts)+gd(ts-1)),'inf');
    err2(j) = norm(u2-(g(ts)+g(ts-1)),'inf');
end
dt = T./Ns;
loglog(dt,err1,dt,err2,'--')
xlabel('\Delta t'); legend('K_1 error','K_2 error','Location','northwest')
hold on;loglog(dt,10*dt.^5,'-.',dt,dt.^3,'-.'); hold off