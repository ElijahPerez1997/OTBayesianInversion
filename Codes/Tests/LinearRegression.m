clc,clear,close all
%%
%Creating Data Set
x=linspace(0,1,500);
sig_star=0.1;%True Sigma
theta1_star=0.5;%True Theta 1
theta2_star=1;%True Theta 2
y=zeros(1,length(x));%Memory prellocation
burn=10000;
for i=1:length(x)
    eps=normrnd(0,sig_star);
    y(i)=theta1_star+theta2_star*x(i)+eps;
end


%%
%MH Within G
n=500;%Number of data points
M=1000000;%Number of samples drawn
a0=1;b0=0.1;%shape and rate of prior
theta1=zeros(1,M);theta1_0=2;theta1(1)=theta1_0;%Initializing theta 1
theta2=zeros(1,M);theta2_0=2;theta2(1)=theta2_0;%Initializing theta 2
theta1_accept=0;%Memory prellocation
theta2_accept=0;%Memory prellocation
a=a0+(n/2);
i=1;j=1;
count=0;%Number of rejected values
r=0.01;%Support of proposal
s=zeros(1,M-1);%Memory prellocation

while i<M
    b=b0+0.5*sum(abs(y-theta1(i)-theta2(i)*x).^2);
    s(i)=gamrnd(a,b);
    
    theta1_s=unifrnd(theta1(i)-(r/2),theta1(i)+(r/2));
    theta2_s=unifrnd(theta2(i)-(r/2),theta2(i)+(r/2));
    
    prior_i=-0.5*[abs(theta1_star-theta1(i)),abs(theta2_star-theta2(i))]*[1/sig_star,0;0,1/sig_star]*[abs(theta1_star-theta1(i));abs(theta2_star-theta2(i))];
    prior_s=-0.5*[abs(theta1_star-theta1_s),abs(theta2_star-theta2_s)]*[1/sig_star,0;0,1/sig_star]*[abs(theta1_star-theta1_s);abs(theta2_star-theta2_s)];
    
    %alpha=s(i)*((0.5*sum(abs(y-theta1(i)-theta2(i)*x).^2))-(0.5*sum(abs(y-theta1_s-theta2_s*x).^2)));%Log of alpha
    %alpha=(prior_s*exp(-0.5*s(i)*sum(abs(y-theta1_s-theta2_s*x).^2)))/(prior_i*exp(-0.5*s(i)*sum(abs(y-theta1(i)-theta2(i)*x).^2)));
    alpha=(((prior_s-0.5*s(i)*sum(abs(y-theta1_s-theta2_s*x).^2)))-((prior_i-0.5*s(i)*sum(abs(y-theta1(i)-theta2(i)*x).^2))));%log of alpha
    u=rand(1);%U hat
    
    
    if log(u)<alpha && log(u)<0
        theta1(i+1)=theta1_s;
        theta2(i+1)=theta2_s;
        theta1_accept(j)=theta1_s;
        theta2_accept(j)=theta2_s;
%       plot(j,theta2_accept(j),'ok')
%       hold on
%       plot(j,theta2_accept(j),'o')
%       hold on
        i=i+1;
        j=j+1;
    else
        theta1(i+1)=theta1(i);
        theta2(i+1)=theta2(i);
        i=i+1;
        count=count+1;
    end
      

end
acceptance_rate=j/i;
f1=figure('Name','Histogram of Theta 1');
histogram(theta1_accept(burn:end),30)
title('Histogram of Theta 1')

f2=figure('Name','Histogram of Theta 2');
histogram(theta2_accept(burn:end),30)
title('Histogram of Theta 2')
exp1=sum(theta1_accept(burn:end))/(j-burn-1);%Expected value of theta 1
exp2=sum(theta2_accept(burn:end))/(j-burn-1);%Expected value of theta 2
lin_reg=exp1+exp2*x;%Linear regression line

f3=figure('Name','Linear Regression Line');
plot(x,y,'ro');
title('Linear Regression Line')
hold on
plot(x,lin_reg,'-k')
legend('Data Points','Regression Line')

f4=figure('Name','Trace Plot of Theta 1');
plot(theta1_accept(burn:end))
title('Trace of Theta 1')

f5=figure('Name','Trace Plot of Theta 2');
plot(theta2_accept(burn:end))
title('Trace of Theta 2')
