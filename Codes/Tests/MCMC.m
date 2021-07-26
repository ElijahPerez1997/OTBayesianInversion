clc,clear,close all
%%
%Creating Data Set
n=500;
x=linspace(0,1,n);
sig_star=0.1;%True Sigma
theta1_star=0.5;%True Theta 1
theta2_star=1;%True Theta 2
y=zeros(1,n);%Memory prellocation

for i=1:n
    eps=normrnd(0,sig_star);
    y(i)=theta1_star+theta2_star*x(i)+eps;
end


%%
%MH Within G

M=1.0E5;%Number of samples drawn
burn=5E4;%Number of values burned
thin=5;
a0=1;b0=0.1;%shape and rate of prior
theta1=zeros(1,M);theta1_0=2;theta1(1)=theta1_0;%Initializing theta 1
theta2=zeros(1,M);theta2_0=2;theta2(1)=theta2_0;%Initializing theta 2
theta1_accept=0;%Memory prellocation
theta2_accept=0;%Memory prellocation
a_aprox=a0+(n/2);
b_aprox=0;
i=1;j=1;
count=0;%Number of rejected values
r=0.01;%Support of proposal
s=zeros(1,M-1);%Memory prellocation

while i<M
    b_aprox=b0+0.5*sum(abs(y-theta1(i)-theta2(i).*x).^2);
    s(i)=gamrnd(a_aprox,1/b_aprox);
    
    theta1_s=unifrnd(theta1(i)-(r/2),theta1(i)+(r/2));
    theta2_s=unifrnd(theta2(i)-(r/2),theta2(i)+(r/2));
    
    prior_i=-0.5*[abs(theta1_star-theta1(i)),abs(theta2_star-theta2(i))]*[1/(sig_star),0;0,1/(sig_star)]*[abs(theta1_star-theta1(i));abs(theta2_star-theta2(i))];
    prior_s=-0.5*[abs(theta1_star-theta1_s),abs(theta2_star-theta2_s)]*[1/(sig_star),0;0,1/(sig_star)]*[abs(theta1_star-theta1_s);abs(theta2_star-theta2_s)];
    
    %alpha=s(i)*((0.5*sum(abs(y-theta1(i)-theta2(i)*x).^2))-(0.5*sum(abs(y-theta1_s-theta2_s*x).^2)));%Log of alpha
    %alpha=(prior_s*exp(-0.5*s(i)*sum(abs(y-theta1_s-theta2_s*x).^2)))/(prior_i*exp(-0.5*s(i)*sum(abs(y-theta1(i)-theta2(i)*x).^2)));
    alpha=(((prior_s-0.5*s(i)*sum(abs(y-theta1_s-theta2_s*x).^2)))-((prior_i-0.5*s(i)*sum(abs(y-theta1(i)-theta2(i)*x).^2))));%log of alpha
    u=rand(1);%U hat
    
    
    if log(u)<alpha && log(u)<0
        theta1(i+1)=theta1_s;
        theta2(i+1)=theta2_s;
        theta1_accept(j)=theta1_s;
        theta2_accept(j)=theta2_s;

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

theta1_burn=theta1(burn:thin:end);
theta2_burn=theta2(burn:thin:end);
s_burn=s(burn:thin:end);
%%

%Joint distribution for Theta 1&2
r=2;
x_star=linspace(0,1,n);
sig_star=0.1;%True Sigma
theta1_star=0.5;%True Theta 1
theta2_star=1;%True Theta 2


V0=0.1*eye(2);
const=1/(2*pi*(det(V0)^0.5));
a0=1;
b0=0.1;
%this is the pdf of the distribution with mean equal to the true value and
%a covariance matrix of 0.1 times the identity
pdf=@(x) const*exp(-0.5*[abs(theta1_star-x(1)),abs(theta2_star-x(2))]*V0*[abs(theta1_star-x(1));abs(theta2_star-x(2))]);

x1=-4:0.1:4;
x2=0.01:0.01:4;
p=zeros(length(x1));
z=zeros(1,length(x2));
for i=1:length(x1)
    for j=1:length(x1)
        p(i,j)=pdf([x1(i);x1(j)]);
    end
end

for i=1:length(x2)
    z(i)=((b0^a0)/gamma(a0))*((1/x2(i))^(a0+1))*exp(-b0/x2(i));
end
f1=figure('Name','Joint prior of theta 1 and theta 2');
mesh(x1,x1,p)
xlabel('Theta 2')
ylabel('Theta 1')
title('Joint Prior of Theta 1 and Theta 2')


%%

%Initializing/Defining Constants
x2d=[ones(length(x_star),1),x_star'];
V=inv(inv(V0)+x2d'*x2d);
mu0=[theta1_star;theta2_star];
mu=V*(inv(V0)*mu0+x2d'*y');
a_star=a0+n/2;
b_star=b0+0.5*(mu0'*inv(V0)*mu0+y*y'-mu'*inv(V)*mu);

%%

%Sigma plots
sig=linspace(0,10*sqrt(b0/(a0+1)),n);

s_exact=1./(sig.^2);
invgamma = @(x,alpha,beta)(gampdf(1./x,alpha,1/beta)./(x.^2));
fh3=figure('Name','Histogram of s');
h=histogram(1./sqrt(s_burn),30);%Histogram of s
title('Histogram of Sigma')
f3=figure('Name','Prior vs Posterior for sigma');
plot(sig,2.*sig.*invgamma(sig.^2,a_star,b_star),'-k');
%xlim([0,1])
hold on;
plot(sig,2.*sig.*invgamma(sig.^2,a0,b0),'-r');
val=h.Values;%Values of the histogram in each bin
bin=zeros(1,length(h.BinEdges)-1);%Memory preallocation
for i=1:length(h.BinEdges)-1
    bin(i)=0.5*(h.BinEdges(i)+h.BinEdges(i+1));%Picking the center value of the bins
end
q=(bin(end)-bin(1))/(length(bin)-1);%Correcting factor for the bin length
plot(bin,val./(sum(val)*q),'b-')
%plot(bin,val,'b-')
hold off;
legend('Posterior','Prior','MCMC');
xlabel('\sigma'); ylabel('marginal density');

%%

%Theta 1 Plots
sx = sqrt(V(1,1)*b_star/a_star);
xx=linspace(0,1,500);
fh1=figure('Name','Histogram of theta 1');
h=histogram(theta1_burn,30);%Histogram of Theta 1
title('Histogram of Theta 1')
f4=figure('Name','Marginal density of theta 1');
plot(xx,tpdf((xx-mu(1))/sx,2*a_star)/sx,'-k');
hold on
sx0 = sqrt(V0(1,1)*b0/a0);
plot(xx,tpdf((xx-mu0(1))/sx0,2*a0)/sx0,'-r');
val=h.Values;%Values of the histogram in each bin
bin=zeros(1,length(h.BinEdges)-1);%Memory preallocation
for i=1:length(h.BinEdges)-1
    bin(i)=0.5*(h.BinEdges(i)+h.BinEdges(i+1));%Picking the center value of the bins
end
q=(bin(end)-bin(1))/(length(bin)-1);%Correcting factor for the bin length
plot(bin,val./(sum(val)*q),'b-')
hold off ; axis tight ;
legend('Posterior','Prior','MCMC');
xlabel(['\theta',num2str(1)]); ylabel('marginal density');

%%

%Theta 2 Plots
xx=linspace(0.5,1.5,500);
fh2=figure('Name','Histogram of theta 2');
h=histogram(theta2_burn,30);%Histogram of Theta 2
title('Histogram of Theta 2')
f5=figure('Name','Marginal density of theta 2');
plot(xx,tpdf((xx-mu(2))/sx,2*a_star)/sx,'-k');
hold on
sx0 = sqrt(V0(1,1)*b0/a0);
plot(xx,tpdf((xx-mu0(2))/sx0,2*a0)/sx0,'-r');

val=h.Values;%Values of the histogram in each bin
bin=zeros(1,length(h.BinEdges)-1);%Memory preallocation
for i=1:length(h.BinEdges)-1
    bin(i)=0.5*(h.BinEdges(i)+h.BinEdges(i+1));%Picking the center value of the bins
end

q=(bin(end)-bin(1))/(length(bin)-1);%Correcting factor for the bin length
plot(bin,val./(sum(val)*q),'b-')
hold off ; axis tight ;
legend('Posterior','Prior','MCMC');
xlabel(['\theta',num2str(2)]); ylabel('marginal density');

%%