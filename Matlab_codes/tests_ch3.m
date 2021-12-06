%%
% In this section we perform some simple numerical experiments using the
% codes described earlier. The whole section has been produced using
% the Matlab command |publish| from the file |tests_ch2.m| 
% that can be downloaded from the book's website.
set(0,'defaultfigureposition',[380 320 540 200],...
'defaultaxeslinewidth',0.9,'defaultaxesfontsize',8,...
'defaultlinelinewidth',1.1,'defaultpatchlinewidth',1.1,...
'defaultlinemarkersize',15); format compact
%%
% Let us first consider the convolution weights for some special choices of
% $K(s)$. If $K(s) = s$, then $K(\partial_t)g = \partial_t g$ and indeed
% the weights of the BDF1 based CQ are as expected (note that here $\Delta t = 1$):
N = 5; T = 5; Ks = @(s) s;
dlt = @(zt) 1-zt;
real(weightsCQ(N,T,Ks,dlt))
%%
% Next, we see that for integration seen as a convolution with $k(t) = 1$
%
% $$K(\partial_t^{\Delta t})g = \int_0^t g(\tau) d\tau, \qquad K(s) = \frac1s,$$ 
%
% the weights of the trapezoid based CQ are
Ks = @(s) 1./s;
dlt = @(zt) 2*(1-zt)./(1+zt);
real(weightsCQ(N,T,Ks,dlt))
%%
% These are just the weights of the compound trapezoid rule with one small
% difference: the final coefficient is 1 and not 1/2. The reason for this
% is that by definition CQ assumes that the input function $g(t)$ is zero
% at $t = 0$.

%%
% As the next example let us compute the fractional integral
% 
% $$I^\alpha g(t) = \frac1{\Gamma(\alpha)}\int_0^t (t-\tau)^{\alpha-1} g(\tau) d\tau = K(\partial_t) g,$$
% 
% where $K(s) = s^{-\alpha}$ and $\alpha > 0$.
alpha = .75; Ks = @(s) s.^(-alpha);
N = 20; T = 1; ts = (0:N)*T/N;
u = evalCQ(@(t) t.^2,Ks,N,T);
%%
% The maximum error is
max(abs(u-(2/gamma(alpha+3))*ts.^(alpha+2)))
%%
% For the final example we choose a hyperbolic symbol
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
% $$u(t) = \sum_{j = 0}^\infty g(t-j).$$
%
% We will use the recursive implementation to solve this system.
% Alternatively we could have applied |evalCQ| with $K^{-1}(s)$.
N = 400; T = 4; ts = (0:N)*T/N;
g = @(t) exp(-100*(t-0.5).^2);
Ks = @(s) 1-exp(-s);
u = solveLT(0,N,g,Ks,T/N);
plot(ts,real(u),ts,g(ts)+g(ts-1)+g(ts-2)+g(ts-3),'--')
legend('computed','exact','Location','west')
%%
% The above computation was done using the trapezoid rule which is
% not an $L$-stable scheme. Using BDF2 with the same parameters, gives a
% pronounced damping effect.
u = evalCQ(g,@(s) 1./Ks(s),N ,T, @(zt) (1-zt)+.5*(1-zt).^2);
plot(ts,real(u),ts,g(ts)+g(ts-1)+g(ts-2)+g(ts-3),'--')
legend('computed','exact','Location','west')
%%
% We now apply the truncated trapezoidal rule. Recall that htis scheme 
% was obtained by a numerical search and is $A$-stable to machine precision.
cs =[0.894431935027534   0.683269582897605   0.610333312441433];
dlt = @(zt) (1-zt)+.5*(1-zt).^2+(1/4)*cs(1)*(1-zt).^3+...
    (1/8)*cs(2)*(1-zt).^4+(1/16)*cs(3)*(1-zt).^5;    

u = evalCQ(g,@(s) 1./Ks(s),N ,T, dlt);
plot(ts,real(u),ts,g(ts)+g(ts-1)+g(ts-2)+g(ts-3),'--')
legend('computed','exact','Location','west')
