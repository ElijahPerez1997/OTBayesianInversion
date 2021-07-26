clc,clear,close all
%You need to include this in your main code



%This part computes the matrix Q
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
XR=1;%This is a stand in quantity do keep Q matrix 2 dimensional 
t=linspace(0,10,1000);
x1=XR/XR(end); %this is the normalized vector of receiver locations (we name it x1)
x2=t/t(end);  %this is normalized time vector (we name it x2)
n1=length(x1); %no. of discrete points in x1
n2=length(x2); %no. of discrete points in x2
n=n1*n2; %the total number of points in 2D signals
lambda=50; %this is the regularizing parameter (we keep it fixed)

delta1=0.5;
s=linspace(-3,3,n);
c=zeros(1,n);
d=2; %this is the dimension of the signals


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
xf=1;
f=zeros(length(xf),n2);
g=zeros(length(xf),n2);

%u=@(x,t,x0,a) 0.5*a*((exp(-100*((x-t-x0-0.5).^2))+exp(-100*((x-t-x0).^2))+exp(-100*((x-t-x0+0.5).^2)))+(exp(-100*((x+t-x0-0.5).^2))+exp(-100*((x+t-x0).^2))+exp(-100*((x+t-x0+0.5).^2))));
W=zeros(1,n);
DSD=zeros(1,n);
L=zeros(1,n);
for i=1:n
f=exp(-((t-4)/delta1).^2)-exp(-((t-5)/delta1).^2)+exp(-((t-6)/delta1).^2);
g=exp(-((t-s(i)-4)/delta1).^2)-exp(-((t-s(i)-5)/delta1).^2)+exp(-((t-s(i)-6)/delta1).^2);
% ff=reshape(f,[n,1]);
% gg=reshape(g,[n,1]);
ff=f';
gg=g';
c=min(min(g));
kk=0.001;
c=abs(c)+kk;
flin=ff+c; flin=flin/(sum(flin));
glin=gg+c; glin = glin/(sum(glin));
    
ds_fg=Sinkhorn2D(flin,glin,Qcell,lambda);
ds_ff=Sinkhorn2D(flin,flin,Qcell,lambda);
ds_gg=Sinkhorn2D(glin,glin,Qcell,lambda);
DSD(i)=abs(ds_fg-0.5*(ds_ff+ds_gg));
L(i)=norm(f-g);
W(i)=Wasserstein(flin',glin',t);
end
plot(s,L/L(1));
ylim([0 2]);
xlim([-3.1 3.1]);
hold on
plot(s,DSD/DSD(1));
plot(s,W/W(1));
xlabel('shift')
ylabel('Distance (normalized)')
legend('L2','DSD','Wass')