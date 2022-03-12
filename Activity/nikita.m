gd0 = [2017 03 20 10 28 0]; % SET UP STARTTIME (UTC)

vtau = 100; % specify voice session duration
dtau = 60;  % specify dataN session duration
stau = 400; % specify dataS session duration

R0 = 6378.137;

Vminday = 14.4;
DMBmo = 4;
SGBmo = 4;

br_m = 4800;
br_f = 1e06;

jd0 = jday(gd0(1),gd0(2),gd0(3),gd0(4),gd0(5),gd0(6)); % find starttime Julian Date
deg = pi/180;

a = 6978.14; e = 0; om = 0; in = 97.792*deg; % specify initial orbit parameters
Om = 35*deg; u = 15*deg;

dOm = 11.125*deg; du = 8*deg; ddu = 12*deg;  % specify constellation parameters
netdata = {dOm, du, ddu, vtau, dtau, stau, Vminday, DMBmo, SGBmo, br_m, br_f};

ini0 = [a e om Om in u];

clear T Tb Tm Ts To s msg succ ans vtau dtau stau HR a e om Om in u dOm du ddu br_m br_f Vminday DMBmo SGBmo
