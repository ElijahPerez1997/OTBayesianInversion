clc,clear,close all
x=[1000,10000,50000,100000,500000];
y=[342.68,6.9463,1.0167,0.05451,0.01446];
x2=(10000000./x.^2);
loglog(x,y,'-')
hold on
loglog(x,x2,'-')
xlabel('M')
ylabel('Wasserstein Distance')
legend('MCMC vs True','O(1/M^{2})')