deg = pi/180;
R0 = 6378.14;
mu = 398600.44;
J2 = 1.08263e-3;

in = astro.sunsyn(300);

alt0 = 550;

p0 = alt0+R0;
T0 = 2*pi*sqrt(p0^3/mu);
v0 = sqrt(mu/p0);

RAANdot_SSO = 360/365; %deg/day
RAANdot = -3*pi*J2*R0^2/p0^2*cosd(in)/T0*86400/deg;

di = astro.sunsyn(550)-astro.sunsyn(300);
DV1 = v0*di*deg*1e3;

r1 = R0+300; r2 = R0+550;
dv1 = sqrt(mu/r1) * (sqrt(2*r2/(r1+r2)) - 1) * 1e3;
dv2 = sqrt(mu/r2) * (1 - sqrt(2*r1/(r1+r2))) * 1e3;

DV2 = sqrt(dv2^2 + DV1^2) + dv1;