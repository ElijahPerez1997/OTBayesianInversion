clc,clear,close all
load ./MatFiles3/Ex_3a_synthetic.mat
nt=length(t);
nr=length(xf);
theta1=0:0.01:0.2;
theta2=3:0.2:7;
was=zeros(length(theta1),length(theta2));
u=@(x,t,x0,a) 0.5*a*((exp(-100*((x-t-x0-0.5).^2))+exp(-100*((x-t-x0).^2))+exp(-100*((x-t-x0+0.5).^2)))+(exp(-100*((x+t-x0-0.5).^2))+exp(-100*((x+t-x0).^2))+exp(-100*((x+t-x0+0.5).^2))));
c=min(min(g));
kk=1.01;
c=abs(c)*kk;  %Linear scaling constant
c=0.01;
for i=1:length(theta1)
    for j=1:length(theta2)
        for k=1:nr
            fk=u(xf(k),t,theta1(i),theta2(j));
            min(fk)
            gk=g(k,:);
            gk=(gk+c)/sum(gk+c);
            fk=(fk+c)/sum(fk+c);
            wk=Wasserstein(fk,gk,t);
            was(i,j)=was(i,j)+wk;
        end
    end
    
end
mesh(was)