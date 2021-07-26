clc,clear,close all
n=500;
t=linspace(0,5,n);
xf=-3:3;
a_star=5;
f=zeros(length(xf),n);
g=zeros(length(xf),n);

h=@(x,x0,a) a*(exp(-100*((x-x0-0.5).^2))+exp(-100*((x-x0).^2))+exp(-100*((x-x0+0.5).^2)));
x0=0;
figure

for j=1:length(xf)
    for i=1:n
        eps1=gamrnd(60,1/60);
        eps2=unifrnd(-0.25,0.25);
        f(j,i)=0.5*h(xf(j)-t(i),x0,a_star)+0.5*h(xf(j)+t(i),x0,a_star);
        g(j,i)=eps1*f(j,i)+eps2;
    end
    subplot(length(xf),1,j)
    plot(t,g(j,:))
    ylim([-5,5])
end
