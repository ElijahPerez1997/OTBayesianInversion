clc,clear,close all
%%
%Creating Data Set
load ./MatFiles/Ex_1a_synthetic.mat

lambda=50;
d=2;
x1=xf/xf(end); %this is the normalized vector of receiver locations (we name it x1)
x2=t/t(end);%this is normalized time vector (we name it x2)
n1=length(x1); %no. of discrete points in x1
n2=length(x2); %no. of discrete points in x2
n=n1*n2; %the total number of points in 2D signals
Qcell=cell(1,d); %we save both Q1 and Q2 matrices as a cell
C1=zeros(n1,n1);
for j=1:n1
    C1(j,:)=(x1-x1(j)).^2;
end
Qcell{1}=exp(-lambda*C1);
C2=zeros(n2,n2);
for j=1:n2
    C2(j,:)=(x2-x2(j)).^2;
end
Qcell{2}=exp(-lambda*C2);
kk=1.01;
c=min(min(g));
c=abs(c)*kk;  %Linear scaling constant
gi=reshape(g,[n,1]);
gi=gi+c;
gi=gi/(sum(gi)); 
DS_gg=Sinkhorn2D(gi,gi,Qcell,lambda);

%%
%MH Within G

M=5E4;%Number of samples drawn
burn=round(0.5*M);%Number of values burned
thin=2;%Thinning period
a0=100;b0=3;%shape and rate of prior for s
theta=zeros(1,M);theta_0=3;theta(1)=theta_0;%Initializing theta 1


%Memory prellocation
a_aprox=a0+1; %shape of posterior of s

i=1;%Index for MC samples
j=1;%j=number of accepted values
count=0;%Number of rejected values

s=zeros(1,M-1);%Memory prellocation
s0=70; %Initial guess for s
s(1)=s0;

kw=2000; %Constant multiplying DSD
fii=zeros(n1,n2);
fss=zeros(n1,n2);


%%

%u=@(x,t,x0,a) 0.5*a*((exp(-100*((x-t-x0-0.5).^2))+exp(-100*((x-t-x0).^2))+exp(-100*((x-t-x0+0.5).^2)))+(exp(-100*((x+t-x0-0.5).^2))+exp(-100*((x+t-x0).^2))+exp(-100*((x+t-x0+0.5).^2))));
   for k=1:n1
      fii(k,:)=Forward_solve(xf(k),t,x0,theta(1));
   end
    
    fi=(fii+c);
    fi=reshape(fi,[n,1]);
    fi=fi/sum(fi);
    DS_fg_i=Sinkhorn2D(fi,gi,Qcell,lambda);
    DS_ff_i=Sinkhorn2D(fi,fi,Qcell,lambda);
    DSD_i=abs(DS_fg_i-0.5*(DS_ff_i+DS_gg));
    DSD_i=kw*DSD_i;
    
    b_aprox=b0+DSD_i;%Rate of posterior of s 
    
while i<M
    disp(i)
    s(i+1)=gamrnd(a_aprox,1/b_aprox);
    
    theta_s=normrnd(theta(i),0.007);%proposal distribution for theta
    
  %u=@(x,t,x0,a) 0.5*a*((exp(-100*((x-t-x0-0.5).^2))+exp(-100*((x-t-x0).^2))+exp(-100*((x-t-x0+0.5).^2)))+(exp(-100*((x+t-x0-0.5).^2))+exp(-100*((x+t-x0).^2))+exp(-100*((x+t-x0+0.5).^2))));
  for k=1:n1
     fss(k,:)=Forward_solve(xf(k),t,x0,theta_s);
  end
    fs=fss+c;
    fs=reshape(fs,[n,1]);
    fs=fs/sum(fs);
    DS_fg_s=Sinkhorn2D(fs,gi,Qcell,lambda);
    DS_ff_s=Sinkhorn2D(fs,fs,Qcell,lambda);
    DSD_s=abs(DS_fg_s-0.5*(DS_ff_s+DS_gg));
    DSD_s=kw*DSD_s;
   
    
    li=-s(i)*kw*DSD_i; %Log liklihood with fi 
    ls=-s(i)*kw*DSD_s; %Log likelihood with fs
    %Note: Contants for the log liklihoods cancle out and therefore are not
    %included
   
  
    alpha=ls-li;%log of alpha
    
    if 2>theta(i) || theta(i)>8
       
       alpha=-1E10;
    end
    if 2>theta_s || theta_s>8
     
       alpha=-1E10;
    end
    
    u=rand(1);%U hat
    
    if log(u)<alpha
        theta(i+1)=theta_s;
      
       %u=@(x,t,x0,a) 0.5*a*((exp(-100*((x-t-x0-0.5).^2))+exp(-100*((x-t-x0).^2))+exp(-100*((x-t-x0+0.5).^2)))+(exp(-100*((x+t-x0-0.5).^2))+exp(-100*((x+t-x0).^2))+exp(-100*((x+t-x0+0.5).^2))));
        
        
        DSD_i=DSD_s;
    
        b_aprox=b0+DSD_i;%Rate of posterior of s  
        i=i+1;
        j=j+1;
    else
        theta(i+1)=theta(i);
        
        i=i+1;
        count=count+1;
    end
      
end
acceptance_rate=j/i;
theta_burn=theta(burn:thin:end);


s_burn=s(burn:thin:end);
%%
nbins=20; %Number of bins in the histogram
figure
plot(theta_burn)
title('Theta Trace')

figure
histogram(theta_burn,nbins)
title('Theta Histogram')