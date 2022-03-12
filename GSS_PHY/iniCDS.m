clear; clc; % Очистка памяти и экрана

h = 60; % Шаг интегрирования
t_fin = 86400 * 90; % Конечное время интегрирования
options = odeset('MaxStep', 90);
dayperiod = 0.1;
global SATDATA PREC STAGESTR

PREC = 2;
STAGESTR = 'N/a\n';
STARTTIME = [2023 1 1 0 0 0];

jd0 = astro.jday(STARTTIME(1), STARTTIME(2), STARTTIME(3), STARTTIME(4), STARTTIME(5), STARTTIME(6));

thrust = 0.005; %NEWTON
Isp = 2500; %SECONDS

XSI = zeros(360,1); %1.496
flag = zeros(360,1);

Sbody = 0.03; %m^2
Sbatt = 0; %m^2
solarmode = 0;

SATDATA = {thrust, Isp, XSI, flag, Sbody, Sbatt, solarmode};

deg = pi/180;

ALT = 415; %km

%      a - e - om - Om - i - u
%orb = [6378.14+ALT 0 0 72.06*deg 97.640225*deg pi/4 26];
orb = [6378.14+ALT 6.225e-4 0 90*deg 51.6417*deg 0 5];
Y = math.orbitize(orb);

clear thrust Isp XSI flag Sm Se solarmode