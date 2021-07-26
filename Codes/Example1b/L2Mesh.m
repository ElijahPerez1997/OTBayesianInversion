clc,clear,close all
load ./MatFiles2/Ex_2a_synthetic.mat
nt=length(t);
nr=length(xf);
m=20;
L2=zeros(m,m);
theta1=linspace(-1,1,m);
theta2=linspace(3,8,m);
c=0;
for i=1:m
    for j=1:m
        u=@(x,t,x0,a) 0.5*a*((exp(-100*((x-t-x0-0.5).^2))+exp(-100*((x-t-x0).^2))+exp(-100*((x-t-x0+0.5).^2)))+(exp(-100*((x+t-x0-0.5).^2))+exp(-100*((x+t-x0).^2))+exp(-100*((x+t-x0+0.5).^2))));
   for k=1:nr
      fi(k,:)=u(xf(k),t,theta1(i),theta2(j));
   end
    
  fi=(fi+c);
 
  Wi_sum=0;
    for l=1:nr
       fi(l,:)=fi(l,:)/sum(fi(l,:)); 
       Wi_sum=Wi_sum+norm(fi(l,:)-g(l,:));
    end
        L2(i,j)=Wi_sum;
       
    end
end

surf(theta1,theta2,L2)
title('L2 Likelihood Mesh')
shading interp
minimum=min(min(L2));
[x,y]=find(L2==minimum);
theta1_s=theta1(x);
theta2_s=theta2(y);