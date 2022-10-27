%%
set(0,'defaultfigureposition',[380 320 540 200],...
'defaultaxeslinewidth',0.9,'defaultaxesfontsize',8,...
'defaultlinelinewidth',1.5,'defaultpatchlinewidth',1.1,...
'defaultlinemarkersize',15); format compact

x = sym('x','real');
Tbdf2 = taylor(1-exp(-1i*x)+.5*(1-exp(-1i*x))^2);
Ttr = taylor(2*(1-exp(-1i*x))/(1+exp(-1i*x)));

%cs =[0.894431935027534   0.683269582897605   0.610333312441433];
%cs = [0.883429591874175   0.653673139874315   0.872603105967837];
%cs = [0.889960187436936   0.692213146163884   0.599920279130717];
cs = [0.893817850529318   0.684154908023834   0.629642997466429];

Tttr = taylor(1-exp(-1i*x)+.5*(1-exp(-1i*x))^2+.25*cs(1)*(1-exp(-1i*x)).^3+(1/8)*cs(2)*(1-exp(-1i*x)).^4+...
                    +(1/16)*cs(3)*(1-exp(-1i*x)).^5);


 dlt_bdf2 = @(zt) 1-zt+.5*(1-zt).^2;
 dlt_tr = @(zt) 2*(1-zt)./(1+zt);
 dlt_ttr = @(zt) 1-zt+.5*(1-zt).^2+(1/4)*cs(1)*(1-zt).^3+(1/8)*cs(2)*(1-zt).^4+...
     (1/16)*cs(3)*(1-zt).^5;
 
 zs = exp(-1i*linspace(0,2*pi,10000));
 plot(dlt_bdf2(zs)); axis equal; hold on
 plot(1e-12+5*1i*linspace(-1,1),'--');
 plot(dlt_ttr(zs),'-.');  hold off
 h = legend('BDF2','TR','TTR');
 set(h,'Interpreter','Latex')
 