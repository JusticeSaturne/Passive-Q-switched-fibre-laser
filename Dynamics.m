clear all
close all
clc
tic
% This script plots differents values of the out put of the function
% PassiveQswitch
Pump1 = 2.5; % Pump power in watts
Pump2 = 10;
[Psp1,time1,N21,Nsa21] = PassiveQswitch(Pump1);
[Psp2,time2,N22,Nsa22] = PassiveQswitch(Pump2);
% The plottings
figure(1)
plot(time1,Psp1,'Linewidth',2);
hold on
plot(time1,Psp2,'r','Linewidth',2);
figure(2)
plot(time1,N21,'Linewidth',2);
hold on
plot(time1,N22,'r','Linewidth',2);
figure(3)
plot(time1,Nsa21,'Linewidth',2);
hold on
plot(time1,Nsa22,'r','Linewidth',2);
toc
