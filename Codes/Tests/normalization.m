clc,clear,close all
n=1000;
t=linspace(0,10,n);
delta1=0.5;
delta2=0.05;
s=-3;
f1=exp(-((t-4)/delta1).^2)-exp(-((t-5)/delta1).^2)+exp(-((t-6)/delta1).^2);
g1=exp(-((t-s-4)/delta1).^2)-exp(-((t-s-5)/delta1).^2)+exp(-((t-s-6)/delta1).^2);
f2=exp(-((t-4)/delta2).^2)-exp(-((t-5)/delta2).^2)+exp(-((t-6)/delta2).^2);
g2=exp(-((t-s-4)/delta2).^2)-exp(-((t-s-5)/delta2).^2)+exp(-((t-s-6)/delta2).^2);
figure
plot(t,f1)
hold on
plot(t,g1)
plot(t,f2)
plot(t,g2)
title('f,g plots with s=-3')
xlabel('time')
legend('f, \delta=0.5','g, \delta=0.5','f, \delta=0.05','g, \delta=0.05')

% WL=zeros(1,n);
% L=WL;
% WE=L;
% WS=L;
% WA=L;
% WLE=L;
% s=linspace(-3,3,n);
% c=zeros(1,n);
% 
% for j=1:2
% for i=1:n
%     
% f1=exp(-((t-4)/delta1).^2)-exp(-((t-5)/delta1).^2)+exp(-((t-6)/delta1).^2);
% g1=exp(-((t-s(i)-4)/delta1).^2)-exp(-((t-s(i)-5)/delta1).^2)+exp(-((t-s(i)-6)/delta1).^2);
% f2=exp(-((t-4)/delta2).^2)-exp(-((t-5)/delta2).^2)+exp(-((t-6)/delta2).^2);
% g2=exp(-((t-s(i)-4)/delta2).^2)-exp(-((t-s(i)-5)/delta2).^2)+exp(-((t-s(i)-6)/delta2).^2);
% if j==1
% f2=f1;
% g2=g1;
% end
%     
%     c(i)=min([f2 g2]);
%     if c(i)<0
%         c(i)=-c(i)+.1;
%         fe=exp(c(i)*f2);
%         ge=exp(c(i)*g2);
%         
%         fle=fe;
%         gle=ge;
%         
%         fa=abs(f2);
%         ga=abs(g2);
%         
%         fs=f2.^2;
%         gs=g2.^2;
%         
%         fL=f2+c(i);
%         gL=g2+c(i);
%     else
%         c(i)=c(i)+.1;
%         fe=exp(c(i)*f2);
%         ge=exp(c(i)*g2);
%         
%         fa=abs(f2);
%         ga=abs(g2);
%         
%         fs=f2.^2;
%         gs=g2.^2;
%         
%         fL=f2+c(i);
%         gL=g2+c(i);
%         
%         fle=fL;
%         gle=gL;
%     end
%     
%     fe=fe/sum(fe);
%     ge=ge/sum(ge);
%     
%     fle=fle/sum(fle);
%     gle=gle/sum(gle);
%     
%     fa=fa/sum(fa);
%     ga=ga/sum(ga);
%     
%     fL=fL/sum(fL);
%     gL=gL/sum(gL);
%     
%     fs=fs/sum(fs);
%     gs=gs/sum(gs);
%     
%     WL(i)=Wasserstein(fL,gL,t);
%     WE(i)=Wasserstein(fe,ge,t);
%     WS(i)=Wasserstein(fs,gs,t);
%     WA(i)=Wasserstein(fa,ga,t);
%     WLE(i)=Wasserstein(fle,gle,t);
%     L(i)=norm(f2-g2);
% end
% 
% 
% figure
% plot(s,L/L(1))
% ylim([0 2]);
% xlim([-3.1 3.1]);
% hold on 
% plot(s,WL/WL(1))
% plot(s,WE/WE(1))
% plot(s,WS/WS(1))
% plot(s,WA/WA(1))
% plot(s,WLE/WLE(1),'r')
% legend('L2','Linear','EXP','Square','Absolute Value','Linear / EXP')
% end