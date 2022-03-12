function [ t1, t2, t3, eph1, eph2, eph3 ] = calcphasing( Da, Di, Yn, deltaRAAN )
%PHASING STAGE calculation 
%   active phasing calculator
%
%   by A.Kharlan, © Yaliny, 2015
%
%     INPUTS:
%     Yn        : nominal orbit parameters for orbit 1
%     Da        : phasing semi-major axis difference, km
%     Di        : phasing inclination difference, degrees
%     deltaRAAN : RAAN difference between orbits, rad
% 
%     OUTPUTS:
%     t1,eph1         : true ephemeris for satellite 1
%     t2,eph2         : true ephemeris for satellite 2
%     t3,eph3         : true ephemeris for satellite 3
% 
%     GLOBAL VARIABLES USED:
%     PROPXSI : yaw XSI-angle defined here ports to the motion model
%     PREC    : defined elsewhere, this determines the modelling precision
%     BOUND     : defined elsewhere, this is sin(ulim) to limit active
%               stage during each orbit. Negative value means phasing,
%               positive value means stationkeeping;

mu = 398600.44;
deg = pi/180;
J2 = 2.634E+10;

h = 25;
options = odeset('MaxStep', 240);

global PROPXSI PREC STAGESTR BOUND;

PREC = 0; %PRECISION FOR ACTIVE STAGE SOLUTION
[ XSIact, deltaOmega, tact, yact, yfilact ] = getactive( Yn, Da, Di );

Dtact = tact(length(tact));

% nom for 1st orb at t0
pn = Yn(1);
l1n = Yn(2);
l2n = Yn(3);
Omn = Yn(4);
inn = Yn(5);
un = Yn(6);
mn = Yn(7);


en = sqrt(l1n^2 + l2n^2);
an = pn * (1-en^2);

driftop = 360/365.25 * deg/86400; % rad/s
Tph = 2*pi*sqrt((an-Da)^3/mu);
driftph = -2*pi*J2/mu/(an-Da)^2*cos(inn-Di*deg)/Tph;

Omk1 = Omn + 360/365.25 * deg/86400 * Dtact;
Om00 = Omk1 - deltaOmega;
Om01 = Omn;
Om02 = Om01 - deltaRAAN;
Om03 = Om01 - 2 * deltaRAAN;

[ Omph2, Omk2, t02, tk2 ] = RAANeq(Om00,0,Om02,driftph,driftop,Dtact,deltaOmega);
[ Omph3, Omk3, t03, tk3 ] = RAANeq(Om00,0,Om03,driftph,driftop,Dtact,deltaOmega);

t_FIN = tk3;

% t02/86400;
% tk2/86400;
% t03/86400;
% tk3/86400;


PREC = 1; %PRECISION FOR EPHEMERIS GENERATION
Y01 = Yn;
    Y01(1) = Y01(1) - Da;
    Y01(5) = Y01(5) - Di*deg;
    Y01(4) = Om00;
Y02 = Yn;
    Y02(1) = Y02(1) - Da;
    Y02(5) = Y02(5) - Di*deg;
    Y02(4) = Omph2;
Y03 = Yn;
    Y03(1) = Y03(1) - Da;
    Y03(5) = Y03(5) - Di*deg;
    Y03(4) = Omph3;
    

fprintf('GENERATING EPHEMERIS FOR ORBIT 1 - active...\n');
PROPXSI = XSIact;
[at1, aeph1] = ode45('earth353_FULL', [0:h:Dtact], Y01, options);
STAGESTR = 'GENERATING EPHEMERIS FOR ORBIT 1 - passive...\n';
[pt1, peph1, peph1fil] = PEG(Dtact, h, t_FIN,0.5,aeph1(length(aeph1),:),Yn,'A');
t1 = [at1;pt1];
eph1 = [aeph1;peph1];


STAGESTR = 'GENERATING EPHEMERIS FOR ORBIT 2 - passive 1...\n';
[pt21, peph21, peph21fil] = PEG(0,h,t02,0.5,Y01,Y01,'A');
fprintf('GENERATING EPHEMERIS FOR ORBIT 2 - active...\n');
PROPXSI = XSIact;
BOUND = -0.819152;
[at2, aeph2] = ode45('earth353_FULL', [t02:h:t02+Dtact], peph21(length(peph21),:), options);
STAGESTR = 'GENERATING EPHEMERIS FOR ORBIT 2 - passive 2...\n';
[pt22, peph22, peph22fil] = PEG(t02+Dtact, h, t_FIN,0.5,aeph2(length(aeph2),:),Yn,'A');
t2 = [pt21;at2;pt22];
eph2 = [peph21; aeph2; peph22];


STAGESTR = 'GENERATING EPHEMERIS FOR ORBIT 3 - passive...\n';
[pt3, peph3, peph3fil] = PEG(0, h, t03,0.5,Y01,Y01,'A');
fprintf('GENERATING EPHEMERIS FOR ORBIT 3 - active...\n');
PROPXSI = XSIact;
BOUND = -0.819152;
[at3, aeph3] = ode45('earth353_FULL', [t03:h:t_FIN], peph3(length(peph3),:), options);
t3 = [pt3;at3];
eph3 = [peph3; aeph3];





end

