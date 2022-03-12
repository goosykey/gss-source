function [A,Pout,Pdis,Pcons,mass,Radiator,Pout_constraint,Pdis_constraint] = RunSpaceY(subar,numelem, MbMo, br)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

clc;

load spaceY_1.mat BIGASSMASK jd map mapdata netdata;
deg = pi/180;
R0 = 6378.14;
au1 = 149597870.700; %km

netdata{8} = MbMo; % MB/month <-
netdata{10} = br; % mbit/s

if ~exist('BIGASSMASK', 'var')
    error('lol wtf u mad bro?');
end

fprintf('Load complete.\n');

%% EPHE

%jd0 = jday(2018,3,20,16,15,0);
EPHE = J2circ([R0+1150 0 53*deg 0],10800,30);
[RR,~] = randv(EPHE(:,2),EPHE(:,3),EPHE(:,6),EPHE(:,4),EPHE(:,3),EPHE(:,7));
Om = EPHE(:,4); in = EPHE(:,3); u = EPHE(:,7);

POI = [33,-84]; % ATLANTA GA 

fprintf('Ephemeris calculation complete.\n');

%% IMPLEM

[ SESDATA, SESCOORD, SESBEG ] = gettrafficearth( jd, 10800, map, mapdata, BIGASSMASK, netdata );

fprintf('Gettraffic done.\n');

[ diameter, sqlength, gain, masscoef ] = spaceY_antenna( subar, numelem );
[ Pout, Nsub] = getener_spacey( 1150, gain, 30, 33, netdata{10}, 60, SESDATA, SESCOORD, SESBEG, jd+1/48, POI);
[ Pdis, Pcons, effcy, mass ] = spaceY_power( Pout, masscoef );

fprintf('---------------------\n\t');
fprintf(char(datetime('now')));
fprintf('   > ANALYSIS RESULTS \n');
fprintf('---------------------\n\n');
fprintf('Round antenna diameter    : %3.2f m\n',diameter);
fprintf('Rect. antenna linear size : %3.2f m\n',sqlength);
fprintf('Antenna basic gain        : %3.2f dB\n',gain);
fprintf('Antenna mass coef.        : %3.4f kg per Watt dissipated\n',masscoef);
fprintf('\n---------------------\n\n');
fprintf('Output RF power   : %3.2f W\n',Pout);
fprintf('Number of subs    : %2.0f\n',Nsub);
fprintf('Dissipated power  : %3.2f W\n',Pdis);
fprintf('Consumed power    : %3.2f W\n',Pcons);
fprintf('Antenna mass      : %3.2f kg\n',mass);
fprintf('\n---------------------\n\n');

% STR = ['good';'bad '];

% CONSTRAINTS
C_Prf = 0.027;
C_Pdis = 1873.64;

Pout_constraint = C_Prf*subar*numelem;
fprintf('RF power constraint   : %8.2f ', Pout_constraint)


if Pout/subar/numelem < C_Prf
    fprintf('...OK\n')
else
    fprintf('...FAILED\n')
end

A = sqlength^2;
Pdis_constraint = C_Pdis*A;
fprintf('Heat constraint       : %8.2f ', Pdis_constraint)

if Pdis/sqlength^2 < C_Pdis
    fprintf('...OK\n')
else
    fprintf('...FAILED\n')
end


% fprintf('RF power constraint : %s\n',STR((Pout/subar/numelem > C_Prf)+1,:));
% fprintf('Heat constraint     : %s\n',STR((Pdis/sqlength^2 > C_Pdis)+1,:));

C_Preal = Pdis/sqlength^2;
Radiator = sqlength^2 * 4 * C_Preal/C_Pdis - sqlength^2;
fprintf('\nAdditional heat dissipation area : %3.2f sq.m\n',Radiator);

fprintf('\n---------------------\n\n');

%% BIG MODELLING

% jt = EPHE(:,1)/86400 + jd;
% gsti = gstime(jt);
% Rsun = sun(jt)*au1;
% 
% global POWERS SUNANGLES;
% 
% nt = numel(jt);
% SUNANGLES = zeros(nt,3);
% POUT = zeros(nt,1);
% PDIS = POUT; PCONS = POUT;
% 
% for i = 1:nt
%     fprintf('|');
%     Recf = R3(gsti(i))*RR(:,i);
%     [ POUT(i), ~] = getener_spacey( 1150, gain, 30, 33, netdata{10}, 60, SESDATA, SESCOORD, SESBEG, jt(i), Recf);
%     [ PDIS(i), PCONS(i), ~, ~ ] = spaceY_power( POUT(i), masscoef );
%     
%     R1sun = Rsun(:,i) - RR(:,i);
%     R1suntsw = xyz2tsw(R1sun, Om(i), in(i), u(i));
%     SUNANGLES(i,1) = acos(R1suntsw(3) / norm(R1suntsw));
%     SUNANGLES(i,2) = atan2(R1suntsw(2),R1suntsw(1));
%     lit = sight(Rsun(:,i),RR(:,i),'e');
%     if strcmp(lit, 'yes')
%         SUNANGLES(i,3)=1;
%     else
%         SUNANGLES(i,3)=0;
%     end
% end
% 
% POWERS = [POUT PDIS PCONS];
% 
% fprintf('\n---------------------\n\n');

end

