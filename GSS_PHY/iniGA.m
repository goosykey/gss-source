clear; clc; % Очистка памяти и экрана

h = 100; % Шаг интегрирования
hL = 0.1;
t_fin = 86400 * 1; % Конечное время интегрирования
options = odeset('MaxStep', 0.1, 'Events','ode_ev');

dayperiod = 0.5;
global SATDATA PREC STAGESTR jdsfin

PREC = 2;
STAGESTR = 'N/a\n';
STARTTIME = [2018 3 20 16 15 0];

jd0 = astro.jday(STARTTIME(1), STARTTIME(2), STARTTIME(3), STARTTIME(4), STARTTIME(5), STARTTIME(6));

jdsfin = jd0*86400 + t_fin;

thrust = 0.00035; %NEWTON
Isp = 2600; %SECONDS

XSI = zeros(360,1); %1.496
flag = zeros(360,1);

Sbody = 0.03; %m^2
Sbatt = 0; %m^2
solarmode = 0;

SATDATA = {thrust, Isp, XSI, flag, Sbody, Sbatt, solarmode};

deg = pi/180;

ALT = 250; %km

%      a - e - om - Om - i - u
orb = [6378.14+ALT 0 0 0*deg 90*deg pi/4 150];
%orb = [6928.14 0 0 10.9*pi/180 94.7924*deg pi/4 650];
Y = math.lambdize(orb);
%Y = math.orbitize(orb);

L50 = Y(6);
Y(6) = jd0*86400;

[L5i, yi] = ode45(@earth775, L50:hL:inf, Y, options);

[ ti,ephi ] = math.delambdize( L5i,yi );

tii = (ti-jd0*86400)/86400;


clear thrust Isp XSI flag Sm Se solarmode L50 L5i yi


