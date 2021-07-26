clc, clear, close all
n=500;
r=2;
x_star=linspace(0,1,n);
sig_star=0.1;%True Sigma
theta1_star=0.5;%True Theta 1
theta2_star=1;%True Theta 2
y=zeros(1,length(x_star));%Memory prellocation
burn=10000;
for i=1:length(x_star)
    eps=normrnd(0,sig_star);
    y(i)=theta1_star+theta2_star*x_star(i)+eps;
end

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
title('Prior of theta 1 and theta 2')
% f2=figure('Name','Prior for sigma');
% plot(x2,z,'-')
% title('Prior for sigma')



x2d=[ones(length(x_star),1),x_star'];
V=inv(inv(V0)+x2d'*x2d);
mu0=[theta1_star;theta2_star];
mu=V*(inv(V0)*mu0+x2d'*y');
a_star=a0+n/2;
b_star=b0+0.5*(mu0'*inv(V0)*mu0+y*y'-mu'*inv(V)*mu);


sig=linspace(0,10*sqrt(b0/(a0+1)),1000);
invgamma = @(x,alpha,beta)(gampdf(1./x,alpha,1/beta)./(x.^2));
f3=figure('Name','Prior vs Posterior for sigma');
plot(sig,2*sig.*invgamma(sig.^2,a_star,b_star),'-k');
xlim([0,1])
hold on;
plot(sig,2*sig.*invgamma(sig.^2,a0,b0),'-r');
hold off;
legend('Posterior','Prior');
xlabel('\sigma'); ylabel('marginal density');


sx = sqrt(V(1,1)*b_star/a_star);
%xx = linspace(mu(1)-4*sx,mu(1)+4*sx,500);
xx=linspace(0,1,500);
f4=figure('Name','Marginal density of theta 1');
plot(xx,tpdf((xx-mu(1))/sx,2*a_star)/sx,'-k');
hold on
sx0 = sqrt(V0(1,1)*b0/a0);
plot(xx,tpdf((xx-mu0(1))/sx0,2*a0)/sx0,'-r');
hold off ; axis tight ;
legend('Posterior','Prior');
xlabel(['\theta',num2str(1)]); ylabel('marginal density');


xx=linspace(0.5,1.5,500);
f5=figure('Name','Marginal density of theta 2');
plot(xx,tpdf((xx-mu(2))/sx,2*a_star)/sx,'-k');
hold on
sx0 = sqrt(V0(1,1)*b0/a0);
plot(xx,tpdf((xx-mu0(2))/sx0,2*a0)/sx0,'-r');
hold off ; axis tight ;
legend('Posterior','Prior');
xlabel(['\theta',num2str(2)]); ylabel('marginal density');