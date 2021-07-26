clc,clear
n=100;
t=linspace(0,5,n);
xf=-3:3;
theta_star=5;
f=zeros(length(xf),n);
g=zeros(length(xf),n);

u=@(x,t,x0,a) 0.5*a*((exp(-100*((x-t-x0-0.5).^2))+exp(-100*((x-t-x0).^2))+exp(-100*((x-t-x0+0.5).^2)))+(exp(-100*((x+t-x0-0.5).^2))+exp(-100*((x+t-x0).^2))+exp(-100*((x+t-x0+0.5).^2))));
x0=0;


for j=1:length(xf)
    
        eps1=gamrnd(60,1/60,1,n);
        eps2=unifrnd(-0.25,0.25,1,n);
        f(j,:)=u(xf(j),t,x0,theta_star);
        g(j,:)=eps1.*f(j,:)+eps2;
    
end
save ./MatFiles/Ex_1a_synthetic.mat g t xf x0 theta_star u
