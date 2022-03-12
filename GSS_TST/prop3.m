clear;
clc;

P = 0.08;
Isp = 1600;
M0 = 1000;
g0 = 9.8066;

R0 = 6378.14;
mu = 398600.44;
r1 = R0 + 600;
r2 = R0 + 35786;

a = (r1 + r2)/2;

dv1 = sqrt(mu/r1) * (sqrt(2*r2/(r1+r2)) - 1);
dv2 = sqrt(mu/r2) * (1 - sqrt(2*r1/(r1+r2)));

dv = (dv1+dv2)*1.15*1e3;

Mp = M0 * (1 - 1/exp(dv/Isp/g0))

mdot = P/Isp/10;
dt_sec = Mp/mdot
dt_hr = dt_sec/3600
dt_day = dt_sec/86400
dt_mo = dt_day/30

vp = sqrt(mu*(2/r1 - 1/a));
va = sqrt(mu*(2/r2 - 1/a));