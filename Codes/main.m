clc,clear,close all
%You need to include this in your main code



%This part computes the matrix Q
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
XR=[-3,-2,-1,0,1,2,3];
t=linspace(0,5,100);
lambda=50; %this is the regularizing parameter (we keep it fixed)
d=2; %this is the dimension of the signals
x1=XR/XR(end); %this is the normalized vector of receiver locations (we name it x1)
x2=t/t(end);  %this is normalized time vector (we name it x2)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Now suppose you have two 2D signals f and g of size n1xn2
%First you need to "correctly" reshape them
%Now suppose you have two signals f and g of size nx1
%Here is how you compute their DSD (dibiased Sinkhorn distance)

%c = ... this is the normalization constant that you need to find
xf=-3:3;
theta_star=[0,5];
x0=theta_star(1);
f=zeros(length(xf),n2);
g=zeros(length(xf),n2);

u=@(x,t,x0,a) 0.5*a*((exp(-100*((x-t-x0-0.5).^2))+exp(-100*((x-t-x0).^2))+exp(-100*((x-t-x0+0.5).^2)))+(exp(-100*((x+t-x0-0.5).^2))+exp(-100*((x+t-x0).^2))+exp(-100*((x+t-x0+0.5).^2))));


%for i=1:100
for j=1:length(xf)
        eps1=normrnd(0,0.1);
        f(j,:)=u(xf(j),t,theta_star(1),theta_star(2));
        g(j,:)=f(j,:)+eps1;
end
ff=reshape(f',[n,1]);
gg=reshape(g',[n,1]);
c=min(min(g));
kk=1.01;
c=abs(c)*kk;
flin=ff+c; flin=flin/(sum(flin));
glin=gg+c; glin = glin/(sum(glin));
    
ds_fg=Sinkhorn2D(flin,glin,Qcell,lambda);
ds_ff=Sinkhorn2D(flin,flin,Qcell,lambda);
ds_gg=Sinkhorn2D(glin,glin,Qcell,lambda);
DSD=abs(ds_fg-0.5*(ds_ff+ds_gg));
plot(1,DSD,'o')
%hold on
%

