clear; clc; % Очистка памяти и экрана

h = 10; % Шаг интегрирования
t_fin = 86400; % Конечное время интегрирования
options = odeset('MaxStep', 90);
dayperiod = 0.5;
global SATDATA PREC STAGESTR

PREC = 2;
STAGESTR = 'N/a\n';
STARTTIME = [2018 1 1 0 0 0];

jd0 = astro.jday(STARTTIME(1), STARTTIME(2), STARTTIME(3), STARTTIME(4), STARTTIME(5), STARTTIME(6));

thrust = 1.2e-3; %NEWTON
Isp = 1200; %SECONDS

XSI = zeros(360,1); %1.496
flag = zeros(360,1);

Sbody = 1; %m^2
Sbatt = 1; %m^2
solarmode = 2;

SATDATA = {thrust, Isp, XSI, flag, Sbody, Sbatt, solarmode};

deg = pi/180;

ALT = 560.993907; %km

%      a - e - om - Om - i - u
orb = [6378.14+ALT 0 0 270.35*deg 97.640225*deg pi/4 200];
%orb = [6928.14 0 0 10.9*pi/180 94.7924*pi/180 pi/4 650];
Y = math.orbitize(orb);

clear thrust Isp XSI flag Sm Se solarmode