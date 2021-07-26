clc,clear,close all
load ./MatFiles3/Ex_3a_synthetic.mat
nt=length(t);
nr=length(xf);
for i=1:nr
   subplot(nr,1,i)
   plot(t,g(i,:))
   axis([0 5 -5 5])
end
