%%
set(0,'defaultfigureposition',[380 320 540 200],...
'defaultaxeslinewidth',0.9,'defaultaxesfontsize',8,...
'defaultlinelinewidth',1.1,'defaultpatchlinewidth',1.1,...
'defaultlinemarkersize',1); format compact
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

A = RKdata(1);
m = size(A,1);
bone = ones(m,1); em = zeros(m,1); em(end) = 1;
dlt= @(zt) inv(A)-zt*inv(A)*(bone*em');

M = 1000;
z = exp(1i*2*pi*(0:M)/(M+1));

egs_2 = [];
for j = 1:M
    egs_2 = [egs_2; eig(dlt(z(j)))];
end
I = imag(egs_2) < 0;
egs_21 = egs_2(I);
egs_21 = sort(egs_21,'ComparisonMethod','real');
egs_21 = [egs_21; conj(flipud(egs_21))];

plot(egs_21,'-'); axis equal;shg; hold on

A = RKdata(0);
m = size(A,1);
bone = ones(m,1); em = zeros(m,1); em(end) = 1;
dlt= @(zt) inv(A)-zt*inv(A)*(bone*em');

M = 1000;
z = exp(1i*2*pi*(0:M)/(M+1));

egs_2 = [];
for j = 1:M
    egs_2 = [egs_2; eig(dlt(z(j)))];
end

I = imag(egs_2) < 0;
egs_21 = egs_2(I);
egs_21 = sort(egs_21,'ComparisonMethod','real');
egs_21 = [egs_21; conj(flipud(egs_21))];
plot(egs_21,'--'); axis equal;shg; 

dlt = @(zt) 1-zt+.5*(1-zt).^2;
plot(dlt(z),'-.');hold off

h = legend('3-stage Radau IIA', '2-stage Radau IIA','BDF2');
set(h,'Interpreter','Latex')

