clear;clc;

P = 2000;%
Isp = 350;

mdot = P/10/Isp;

dV_re = Isp*10*log(100/50);

dt = 50/mdot;

V_bo = dV_re - 10*dt
