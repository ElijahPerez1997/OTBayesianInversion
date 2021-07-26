clc,clear,close all
load ./MatFiles2/Ex_2a_synthetic.mat

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
m=20;
DSD=zeros(m,m);
theta1=linspace(-1,1,m);
theta2=linspace(3,8,m);
kw=1;
for i=1:m
    for j=1:m
        for k=1:n1
            fii(k,:)=Forward_solve(xf(k),t,theta1(i),theta2(j));
        end
        fi=(fii+c);
        fi=reshape(fi,[n,1]);
        fi=fi/sum(fi);
        DS_fg_i=Sinkhorn2D(fi,gi,Qcell,lambda);
        DS_ff_i=Sinkhorn2D(fi,fi,Qcell,lambda);
        DSD(i,j)=abs(DS_fg_i-0.5*(DS_ff_i+DS_gg));
    end
end
surf(theta1,theta2,kw*DSD)
title('Sinkhorn Likelihood Mesh')
shading interp
minimum=min(min(DSD));
[x,y]=find(DSD==minimum);
theta1_s=theta1(x);
theta2_s=theta2(y);