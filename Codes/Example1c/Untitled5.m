clc,clear,close all
t=linspace(0,6,100000);
true1=normpdf(t,0.1,0.001);
true2=normpdf(t,5,0.01);
plot(t,true1)
hold on
plot(t,true2)
Wasserstein(true1,true2,t)