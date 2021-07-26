clc,clear,close all
%%
%Creating Data Set
load ./MatFiles/Ex_1a_synthetic.mat
nt=length(t);
nr=length(xf);
%%
%MH Within G

M=9E5;%Number of samples drawn
burn=1E4;%Number of values burned
thin=2;%Thinning period
a0=1;b0=0.10;%shape and rate of prior
theta1=zeros(1,M);theta1_0=3;theta1(1)=theta1_0;%Initializing theta 1

theta1_accept=0;%Memory prellocation
%Memory prellocation
a_aprox=a0+nt; 
%a_aprox=a0+(n/2);
b_aprox=0;
i=1;j=1;
count=0;%Number of rejected values

s=zeros(1,M-1);%Memory prellocation
s0=70;
s(1)=s0;
kk=1.1;
fi=zeros(nr,nt);
fs=zeros(nr,nt);
c=min(min(g));
c=abs(c)*kk;
gi=g+c;
for l=1:nr
    gi(l,:)=gi(l,:)/sum(gi(l,:)); 
end



while i<M
    disp(i)
    u=@(x,t,x0,a) 0.5*a*((exp(-100*((x-t-x0-0.5).^2))+exp(-100*((x-t-x0).^2))+exp(-100*((x-t-x0+0.5).^2)))+(exp(-100*((x+t-x0-0.5).^2))+exp(-100*((x+t-x0).^2))+exp(-100*((x+t-x0+0.5).^2))));
   for k=1:nr
      fi(k,:)=u(xf(k),t,x0,theta1(i));
   end
    
  fi=(fi+c);
 
  Wi_sum=0;
    for l=1:nr
       fi(l,:)=fi(l,:)/sum(fi(l,:)); 
       Wi_sum=Wi_sum+Wasserstein(fi(l,:),gi(l,:),t);
    end
    
    b_aprox=b0+Wi_sum;  
    s(i+1)=gamrnd(a_aprox,1/b_aprox);
    
    theta1_s=normrnd(theta1(i),0.001);%proposal distributions for theta 1
 
  for k=1:nr
     fs(k,:)=u(xf(k),t,x0,theta1_s);
  end
    fs=fs+c;
    
    Ws_sum=0;
    for l=1:nr
       fs(l,:)=fs(l,:)/sum(fs(l,:)); 
       Ws_sum=Ws_sum+Wasserstein(fs(l,:),gi(l,:),t);
    end
    
    li=-s(i)*Wi_sum; %Log liklihood with fi 
    ls=-s(i)*Ws_sum; %Log likelihood with fs
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