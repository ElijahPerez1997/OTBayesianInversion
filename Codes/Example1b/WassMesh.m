clc, clear,close all
load ./MatFiles2/Ex_2a_synthetic.mat
nt=length(t);
nr=length(xf);
c=min(min(g));
kk=1.01;
c=abs(c)*kk;  %Linear scaling constant
gi=g+c;
for l=1:nr
    gi(l,:)=gi(l,:)/sum(gi(l,:)); 
end
m=20;
Wsum=zeros(m,m);
theta1=linspace(-1,1,m);
theta2=linspace(3,8,m);

for i=1:m
    for j=1:m
            u=@(x,t,x0,a) 0.5*a*((exp(-100*((x-t-x0-0.5).^2))+exp(-100*((x-t-x0).^2))+exp(-100*((x-t-x0+0.5).^2)))+(exp(-100*((x+t-x0-0.5).^2))+exp(-100*((x+t-x0).^2))+exp(-100*((x+t-x0+0.5).^2))));
            for k=1:nr
                fs(k,:)=u(xf(k),t,theta1(i),theta2(j));
            end
            fs=fs+c;
    
            W=0;
            for l=1:nr
                fs(l,:)=fs(l,:)/sum(fs(l,:)); 
                W=W+Wasserstein(fs(l,:),gi(l,:),t);
            end
            Wsum(i,j)=W;
    end
end

surf(theta1,theta2,Wsum)
title('Wasserstein Likelihood Mesh')
shading interp
minimum=min(min(Wsum));
[x,y]=find(Wsum==minimum);
theta1_s=theta1(x);
theta2_s=theta2(y);