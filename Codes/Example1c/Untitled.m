clc,clear,close all
load ./MatFiles3/Ex_3a_synthetic.mat
nt=length(t);
nr=length(xf);
f=zeros(nr,nt);

figure
for j=1:nr
    for i=1:nt
        f(j,i)=g(j,i);
    end
    subplot(length(xf),1,j)
    plot(t,f(j,:))
   ylabel("g(t,"+char(sym("x_"+j))+")")
    ylim([-5,5])
end
xlabel('t')