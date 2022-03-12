clear; clc; % Очистка памяти и экрана

h = 25; % Шаг интегрирования
t_fin = 86400 * 365; % Конечное время интегрирования
options = odeset('MaxStep', 90);
dayperiod = 0.5;
global PROPXSI PREC STAGESTR STARTTIME BOUND SATNAME
PROPXSI = 0;
PREC = 2;
STAGESTR = 'N/a\n';
STARTTIME = [2017 1 1 0 0 0];
BOUND = 0;
SATNAME = 'SL';

deg = pi/180;

ALT = 600; %km

%      a - e - om - Om - i - u
orb = [42164*0.43 0 0 10*deg 5*deg 0 3500];
%orb = [6928.14 0 0 10.9*pi/180 94.7924*pi/180 pi/4 650];
Y = orbitize(orb);
