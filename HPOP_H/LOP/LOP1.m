 
clear; clc;

tic;

STARTTIME = [2023 1 1 0 0 0];

jd0 = astro.jday(STARTTIME(1), STARTTIME(2), STARTTIME(3), STARTTIME(4), STARTTIME(5), STARTTIME(6));

deg = pi/180;

ALT = 250; %km

%      a - e - om - Om - i - u
%orb = [6378.14+ALT 0 0 72.06*deg 97.640225*deg pi/4 26];
orb = [6378.14+ALT 0 0 0 astro.sunsyn(ALT)*deg pi/4 200];
Y = math.orbitize(orb);

options = odeset('Events', @event_SMA);

[ ti, yi ] = ode45( @(t,y)earth353_LOP(t,y,'jd0', jd0), 0:1e8, Y, options);



toc;

p = yi(:,1);
l1 = yi(:,2);
l2 = yi(:,3);
e = sqrt(l1.^2 + l2.^2);
a = p./(1-e.^2);
plot(ti/86400, a)
 % 
LOP/LOP2.m
 
clear; clc;

tic;

STARTTIME = [2023 1 1 0 0 0];

jd0 = astro.jday(STARTTIME(1), STARTTIME(2), STARTTIME(3), STARTTIME(4), STARTTIME(5), STARTTIME(6));

deg = pi/180;

ALT = 415; %km

%      a - e - om - Om - i - u

orb = [6378.14+ALT 0.0016 0 0 astro.sunsyn(ALT)*deg pi/4 200];

[R0,V0] = math.randv(orb(1),orb(2),orb(5),orb(4),orb(3),orb(6)-orb(3));

%options = odeset('MaxStep', 30,'RelTol', 1e-05);
options = odeset('Events', @event_R, 'MaxStep', 90);

[ ti, yi ] = ode45( @(t,y)earth324_LOP(t,y,'jd0', jd0), 0:1e8, [R0;V0;200], options);



toc;

rabs = math.absVec(yi(:,1:3),2);


plot(ti/86400, rabs-6378.14);
hold on; grid on;

xlabel('time, days');
ylabel('spheroid altitude, km');

%% ANALYTICAL

[am, em, omm, Omm, inm, Mm, rm] = meanElementsP(orb(1), orb(2), orb(3), orb(4), orb(5), orb(6));
orbm = [am, em, omm, Omm, inm, Mm, rm];

tic

options = odeset('Events', @event_R_analyt, 'maxStep', 3600*3);
[ta, ya] = ode45(@(t,y)earthLOP(t,y), 0:86400:86400*365*30, orb(1:6)',options );

toc

plot(ta/86400,ya(:,1)-6378.14, 'linewidth', 2); hold on; grid on
 % 
LOP/LOP_DEMO.m
 
clear; clc;

tic;

STARTTIME = [2023 1 1 0 0 0];

jd0 = astro.jday(STARTTIME(1), STARTTIME(2), STARTTIME(3), STARTTIME(4), STARTTIME(5), STARTTIME(6));

deg = pi/180;

ALT = 415; %km

%      a - e - om - Om - i - u

orb = [6378.14+ALT 0.0016 0 0 astro.sunsyn(ALT)*deg pi/4 200];

%% ANALYTICAL

tic

options = odeset('Events', @event_R_analyt, 'maxStep', 3600*3);
[ta, ya] = ode45(@(t,y)earthLOP(t,y), 0:86400:86400*365*30, orb(1:6)',options );

toc

plot(ta/86400,ya(:,1)-6378.14, 'linewidth', 2); hold on; grid on
