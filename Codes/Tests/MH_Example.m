clc,clear%close all
theta0=1;%Initial value of theta
a=3;%Alpha in beta distribution
b=3;%Beta in beta distribution
sumy=10;%Sum of data points
M=200000;%Number of samples drawn
burn=5000;%Burn in 
n=30;%Number of observations

r=0.1;%Support of proposal
theta=zeros(M,1);%Memory preallocation
theta(1)=theta0;
i=1;j=1;
theta_accept=0;%Memory preallocation
count=0;%Number of rejected samples 


%f1=figure('Name','Markov chain accepted samples');
while i<M
    acceptance=rand(1);%U hat
    
    theta_s=mod(unifrnd(theta(i)-(r/2),theta(i)+(r/2)),1);
    
   
   alpha=(theta_s^(sumy+a-1)*(1-theta_s)^(n-sumy+b-1))/(theta(i)^(sumy+a-1)*(1-theta(i))^(n-sumy+b-1));
   if alpha>=acceptance 
       theta(i+1)=theta_s;
       theta_accept(j)=theta_s;
       %plot(j,theta_accept(j),'ok')
       %hold on
       j=j+1;
       i=i+1;
 
   else
       theta(i+1)=theta(i);
       %plot(i+1,theta(i+1),'*r')
       i=i+1;
       count=count+1;
   end
   
end
% hold off
accpetance_rate=(M-count)/M;
f2=figure('Name','Normalized Histogram vs. True Distribution');
h=histogram(theta_accept(burn:end),30,'Normalization','pdf');
val=h.Values;%Values of the histogram in each bin
bin=zeros(1,length(h.BinEdges)-1);%Memory preallocation
for i=1:length(h.BinEdges)-1
    bin(i)=0.5*(h.BinEdges(i)+h.BinEdges(i+1));%Picking the center value of the bins
end

q=(bin(end)-bin(1))/(length(bin)-1);%Correcting factor for the bin length

plot(bin,val./(sum(val)*q),'r-')

hold on
plot(bin,val,'ok')
x=0:0.01:1;
plot(x,betapdf(x,13,23))