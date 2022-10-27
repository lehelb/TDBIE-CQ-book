%%
% In this section we perform some simple numerical experiments using the
% codes described earlier. The whole section has been produced using
% the Matlab command |publish| from the file |tests_ch8.m| 
% that can be downloaded from the book's website.
set(0,'defaultfigureposition',[380 320 540 200],...
'defaultaxeslinewidth',0.9,'defaultaxesfontsize',10,...
'defaultlinelinewidth',1.1,'defaultpatchlinewidth',1.1,...
'defaultlinemarkersize',15); format compact
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
lw = 'linewidth'
%%
% Let us first consider the symbol $$K(s) = e^{-s^\beta}$$ with $\beta \in [0,1]$. 
% For $\beta = 0$, $K(\partial_t)g = g$, for $\beta = 1$ we get the shift
% $K(\partial_t)g(t) = g(t-1)$. For a particular choice of $g$, we also
% plot the result for different values of $\beta$ between these two
% extreme values.
N = 400; T = 4; ts = (0:N)*T/N;
g = @(t) exp(-400*(t-0.5).^2);
Ks = @(s) exp(-s.^(1/2));
u1 = evalRKCQ(g,Ks,N,T,1); %3-stage Radau IIA method
Ks = @(s) exp(-s.^(0.9));
u2 = evalRKCQ(g,Ks,N,T,1); %3-stage Radau IIA method

plot(ts,exp(-1)*g(ts),'--',ts, u1, '-.',ts, u2, ts, g(ts-1),':',lw,2)
h = legend('$\beta = 0$','$\beta = 1$','$\beta = 1.8$','$\beta = 2$','Location','east');
