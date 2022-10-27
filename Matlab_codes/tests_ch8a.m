%%
set_param_cq;
%%
T = 4; N = 400;
ts = (0:N)*T/N;
beta = 1; 
[us,us_ex,NQ1] = kernel_beta(2,beta,T,N);
figure(1);
semilogy(ts,abs(us-us_ex))

[us,us_ex,NQ2] = kernel_beta(4,beta,T,N);
hold on; semilogy(ts,abs(us-us_ex),'--'); hold off
xlabel('$t$','interpreter','Latex'); ylabel('error','interpreter','Latex')
legend(['$M = $' num2str(NQ1)],['$M = $' num2str(NQ2)])
saveas(gcf,'html/kernel_beta1','epsc')

beta = 1.8; 
[us,us_ex,NQ1] = kernel_beta(4,beta,T,N);
figure(2)
semilogy(ts,abs(us-us_ex))
xlabel('$t$','interpreter','Latex'); ylabel('error','interpreter','Latex')


[us,us_ex,NQ2] = kernel_beta(8,beta,T,N);
hold on; semilogy(ts,abs(us-us_ex),'--');
[us,us_ex,NQ3] = kernel_beta(16,beta,T,N);
semilogy(ts,abs(us-us_ex),'-.');hold off
legend(['$M = $' num2str(NQ1)],['$M = $' num2str(NQ2)],['$M = $' num2str(NQ3)])
saveas(gcf,'html/kernel_beta1p8','epsc')

%%
T = 4; N = 400;
ts = (0:N)*T/N;
r = .3; 
[us,us_ex,NQ1] = kernel_2d(2,r,T,N);
figure(1);
semilogy(ts,abs(us-us_ex)); hold on
[us,us_ex,NQ2] = kernel_2d(4,r,T,N);
semilogy(ts,abs(us-us_ex),'--'); hold off
legend(['$M = $' num2str(NQ1)],['$M = $' num2str(NQ2)])
saveas(gcf,'html/kernel_2d','epsc')
