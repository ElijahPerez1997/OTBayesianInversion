clc,clear,close all
n=100;
t=linspace(0,5,n);
xf=-3:3;
theta_star=[0,5];
x0=theta_star(1);
f=zeros(length(xf),n);
g=zeros(length(xf),n);

u=@(x,t,x0,a) 0.5*a*((exp(-100*((x-t-x0-0.5).^2))+exp(-100*((x-t-x0).^2))+exp(-100*((x-t-x0+0.5).^2)))+(exp(-100*((x+t-x0-0.5).^2))+exp(-100*((x+t-x0).^2))+exp(-100*((x+t-x0+0.5).^2))));



for j=1:length(xf)
    for i=1:length(t)
        eps1=normrnd(0,0.1);
        f(j,i)=u(xf(j),t(i),theta_star(1),theta_star(2));
        g(j,i)=f(j,i)+eps1;
    end
end
save ./MatFiles2/Ex_2a_synthetic.mat g t xf theta_star