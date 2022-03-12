clear; clc; % ������� ������ � ������

h = 50; % ��� ��������������
t_fin = 86400 * 365; % �������� ����� ��������������
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

Sbody = 0.7; %m^2
Sbatt = 1.2; %m^2
solarmode = 2;

SATDATA = {thrust, Isp, XSI, flag, Sbody, Sbatt, solarmode};

deg = pi/180;

ALT = 350; %km

%      a - e - om - Om - i - u
%orb = [6378.14+ALT 0 0 72.06*deg 97.640225*deg pi/4 26];
orb = [6378.14+ALT 0 0 0 astro.sunsyn(ALT)*deg pi/4 200];
Y = math.orbitize(orb);

clear thrust Isp XSI flag Sm Se solarmode