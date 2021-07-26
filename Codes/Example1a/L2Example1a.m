clc,clear,close all
%%
%Creating Data Set
n=100;
t=linspace(0,5,n);
xf=-3:3;
theta_star=5;
f=zeros(length(xf),n);
g=zeros(length(xf),n);

h=@(x,x0,a) a*(exp(-100*((x-x0-0.5).^2))+exp(-100*((x-x0).^2))+exp(-100*((x-x0+0.5).^2)));
x0=0;


for j=1:length(xf)
    for i=1:n
        eps1=gamrnd(60,1/60);
        eps2=unifrnd(-0.25,0.25);
        f(j,i)=0.5*h(xf(j)-t(i),x0,theta_star)+0.5*h(xf(j)+t(i),x0,theta_star);
        g(j,i)=eps1*f(j,i)+eps2;
    end

end

%%
%MH Within G

M=5E4;%Number of samples drawn
burn=1E4;%Number of values burned
thin=4;
a0=1;b0=0.10;%shape and rate of prior
theta1=zeros(1,M);theta1_0=3;theta1(1)=theta1_0;%Initializing theta 1

theta1_accept=0;%Memory prellocation
%Memory prellocation
a_aprox=a0+(n/2); 

b_aprox=0;
i=1;j=1;
count=0;%Number of rejected values

s=zeros(1,M-1);%Memory prellocation
s0=70;
s(1)=s0;
kk=1.1;
fi=zeros(length(xf),n);
fs=zeros(length(xf),n);
while i<M
  for k=1:length(xf)
    for j=1:n
     fi(k,j)=0.5*h(xf(k)-t(j),x0,theta1(i))+0.5*h(xf(k)+t(j),x0,theta1(i));
    end
  end
   
%      c=min(min(g));
%      c=abs(c)*kk;
c=0;
    fi=(fi+c);
    gi=g+c;

    
    fi_sum=0;
    for l=1:length(xf)
       %gi(l,:)=gi(l,:)/sum(gi(l,:));
       %fi(l,:)=fi(l,:)/sum(fi(l,:)); 
       fi_sum=fi_sum+norm(fi(l,:)-gi(l,:));
    end
    
    
    b_aprox=b0+0.5*fi_sum;  
    s(i+1)=gamrnd(a_aprox,1/b_aprox);
    
    
    theta1_s=normrnd(theta1(i),sqrt(0.5));%proposal distributions for theta 1
 
    
  for k=1:length(xf)
    for j=1:n
     fs(k,j)=0.5*h(xf(k)-t(j),x0,theta1_s)+0.5*h(xf(k)+t(j),x0,theta1_s);
    end
  end
    fs=fs+c;
    
    fs_sum=0;
    for l=1:length(xf)
      % gi(l,:)=gi(l,:)/sum(gi(l,:));
       %fs(l,:)=fs(l,:)/sum(fs(l,:)); 
       fs_sum=fs_sum+norm(fs(l,:)-gi(l,:));
    end
    
    
    li=-s(i)*0.5*fi_sum; %Log liklihood with fi 
    ls=-s(i)*0.5*fs_sum; %Log likelihood with fs
    %Note: Contants for the log liklihoods cancle out and therefore are not
    %included
    

    prior_i=1/6;
    prior_s=1/6;
    alpha=(((prior_s+ls))-((prior_i+li)));%log of alpha
    if 2>theta1(i) || theta1(i)>8
       prior_i=0; 
       alpha=-1E10;
    end
    if 2>theta1_s || theta1_s>8
       prior_S=0;
       alpha=-1E10;
    end
    
    
   
    u=rand(1);%U hat
    
    
    if log(u)<alpha && log(u)<0
        theta1(i+1)=theta1_s;
        
        theta1_accept(j)=theta1_s;
        
        
        i=i+1;
        j=j+1;
    else
        theta1(i+1)=theta1(i);
        
        i=i+1;
        count=count+1;
    end
      
end
acceptance_rate=j/i;

theta1_burn=theta1(burn:thin:end);

s_burn=s(burn:thin:end);
%%
nbins=20; %Number of bins in the histogram
figure
plot(theta1_burn)
title('Theta Trace')

figure
histogram(theta1_burn, nbins)
title('Theta Histogram')
