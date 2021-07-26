clc,clear,close all
nt=100;
theta1=linspace(-4,4,nt);

theta1_star=2;
theta2_star=1;
n=500;
x=linspace(0,1,n);
sig_star=0.1;
eps=normrnd(0,sig_star,n,1)';
W=zeros(1,nt);
L=W;
g=theta1_star+theta2_star*x+eps;
for i=1:nt
    f=theta1(i)+theta2_star*x;
    c=min([f g]);
    if c<0
        c=-c+1;
        fi=(f+c);
        gi=g+c;
    else
        fi=f;
        gi=g;
    end
    fi=fi/sum(fi);
    gi=gi/sum(gi);
    W(i)=Wasserstein(fi,gi,x);
    L(i)=norm(f-g);
end
figure
plot(theta1,W,'-o')
figure
plot(theta1,L,'v-')