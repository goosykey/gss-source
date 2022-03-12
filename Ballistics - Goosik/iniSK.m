clear; clc; % Очистка памяти и экрана

h = 25; % Шаг интегрирования
t_fin = 86400 * 365 *  1; % Конечное время интегрирования
options = odeset('MaxStep', 90);
dayperiod = 0.5;
global PROPXSI PREC STAGESTR STARTTIME BOUND SATNAME
PROPXSI = 1.496;
PREC = 2;
STAGESTR = 'N/a\n';
STARTTIME = [2017 3 20 10 28 0];
% BOUND = 0.984808;
BOUND = 0;
SATNAME = 'yasat';

deg = pi/180;

ALT = 500; %km

%      a - e - om - Om - i - u
orb = [6378.14+ALT 0 0 0*deg sunsyn(ALT)*deg pi/4 60];
%orb = [6928.14 0 0 10.9*pi/180 94.7924*pi/180 pi/4 650];
Y = orbitize(orb);