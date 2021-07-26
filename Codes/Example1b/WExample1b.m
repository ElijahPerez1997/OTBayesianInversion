clc,clear,close all
%%
%Creating Data Set
load ./MatFiles2/Ex_2a_synthetic.mat
nt=length(t);
nr=length(xf);
%%
%MH Within G

M=5E4;%Number of samples drawn
burn=round(0.5*M);%Number of values burned
thin=4;%Thinning period
a0=1200;b0=2;%shape and rate of prior for s
theta1=zeros(1,M);theta1_0=0.6;theta1(1)=theta1_0;%Initializing theta 1
theta2=zeros(1,M);theta2_0=3;theta2(1)=theta2_0;%Initializing theta 2

%Memory prellocation
a_aprox=a0+nr; %shape of posterior of s

b_aprox=0;
i=1;j=1;% j=number of accepted values
count=0;%Number of rejected values

s=zeros(1,M-1);%Memory prellocation
s0=70; %Initial guess for s
s(1)=s0;
kk=1.01;
kw=5; %Constant multiplying wasserstein distance
fi=zeros(nr,nt);
fs=zeros(nr,nt);
c=min(min(g));
c=abs(c)*kk;  %Linear scaling constant
gi=g+c;
for l=1:nr
    gi(l,:)=gi(l,:)/sum(gi(l,:)); 
end

u=@(x,t,x0,a) 0.5*a*((exp(-100*((x-t-x0-0.5).^2))+exp(-100*((x-t-x0).^2))+exp(-100*((x-t-x0+0.5).^2)))+(exp(-100*((x+t-x0-0.5).^2))+exp(-100*((x+t-x0).^2))+exp(-100*((x+t-x0+0.5).^2))));
   for k=1:nr
      fi(k,:)=u(xf(k),t,theta1(1),theta2(1));
   end
    
  fi=(fi+c);
 
  Wi_sum=0;
    for l=1:nr
       fi(l,:)=fi(l,:)/sum(fi(l,:)); 
       Wi_sum=Wi_sum+Wasserstein(fi(l,:),gi(l,:),t);
    end
    Wi_sum=kw*Wi_sum;
    b_aprox=b0+Wi_sum;%Rate of posterior of s  


while i<M
    disp(i)
    s(i+1)=gamrnd(a_aprox,1/b_aprox);
    
    thetas=mvnrnd([theta1(i),theta2(i)],[0.00004,0.00004]);%proposal distribution for theta 1
    theta1_s=thetas(:,1);
    theta2_s=thetas(:,2);
    u=@(x,t,x0,a) 0.5*a*((exp(-100*((x-t-x0-0.5).^2))+exp(-100*((x-t-x0).^2))+exp(-100*((x-t-x0+0.5).^2)))+(exp(-100*((x+t-x0-0.5).^2))+exp(-100*((x+t-x0).^2))+exp(-100*((x+t-x0+0.5).^2))));
  for k=1:nr
     fs(k,:)=u(xf(k),t,theta1_s,theta2_s);
  end
    fs=fs+c;
    
    Ws_sum=0;
    for l=1:nr
       fs(l,:)=fs(l,:)/sum(fs(l,:)); 
       Ws_sum=Ws_sum+Wasserstein(fs(l,:),gi(l,:),t);
    end
    Ws_sum=kw*Ws_sum;
    li=-s(i)*kw*Wi_sum; %Log liklihood with fi 
    ls=-s(i)*kw*Ws_sum; %Log likelihood with fs
    %Note: Contants for the log liklihoods cancle out and therefore are not
    %included
   
  
    alpha=ls-li;%log of alpha
    if -3>theta1(i) || theta1(i)>3
        alpha=-1E10;
    end
    if -3>theta1_s || theta1_s>3
        alpha=-1E10;
    end 
    if 2>theta2(i) || theta2(i)>8
       
       alpha=-1E10;
    end
    if 2>theta2_s || theta2_s>8
     
       alpha=-1E10;
    end
    
    u=rand(1);%U hat
    
    if log(u)<alpha && log(u)<0
        theta1(i+1)=theta1_s;
        theta2(i+1)=theta2_s;
%    u=@(x,t,x0,a) 0.5*a*((exp(-100*((x-t-x0-0.5).^2))+exp(-100*((x-t-x0).^2))+exp(-100*((x-t-x0+0.5).^2)))+(exp(-100*((x+t-x0-0.5).^2))+exp(-100*((x+t-x0).^2))+exp(-100*((x+t-x0+0.5).^2))));
%    for k=1:nr
%       fi(k,:)=u(xf(k),t,theta1(i+1),theta2(i+1));
%    end
%     
%   fi=(fi+c);
%  
%   Wi_sum=0;
%     for l=1:nr
%        fi(l,:)=fi(l,:)/sum(fi(l,:)); 
%        Wi_sum=Wi_sum+Wasserstein(fi(l,:),gi(l,:),t);
%     end
%     Wi_sum=kw*Wi_sum;
Wi_sum=Ws_sum;
    b_aprox=b0+Wi_sum;%Rate of posterior of s  
       
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
nbins=20; %Number of bins in each histogram
figure
plot(theta1_burn)
title('Theta 1 Trace')

figure
histogram(theta1_burn,nbins)
title('Theta 1 Histogram')

figure
plot(theta2_burn)
title('Theta 2 Trace')

figure
histogram(theta2_burn,nbins)
title('Theta 2 Histogram')
