clear; clc; % Очистка памяти и экрана

h = 10; % Шаг интегрирования
t_fin = 86400 * 2; % Конечное время интегрирования
options = odeset('MaxStep', 240);
dayperiod = 0.5;
global PROPXSI PREC STAGESTR STARTTIME BOUND SATNAME
PROPXSI = -1;
PREC = 0;
STAGESTR = 'N/a\n';
STARTTIME = [2017 3 20 10 28 0];
BOUND = 0.984808;
SATNAME = 'yasat';

deg = pi/180;


%      a - e - om - Om - i - u
orb = [6978.14 0 0 10.9*pi/180 97.7924*pi/180 pi/4 650];
%orb = [6928.14 0 0 10.9*pi/180 94.7924*pi/180 pi/4 650];
Y = orbitize(orb);
