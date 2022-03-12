function [ STATE510, CART510, LATLON510 ] = seedconstel( PRIMARY_STATE, jd, pla, spp, spreadangle )
%GET ORBITAL POSITIONS OF 510 SATS
%   PLANE #1 is PRIMARY
%   RETURNS CIRCULAT POSITIONS, IGNORES ECCENTRICITY

deg = pi/180;
angleraw = [90,155,174, 176];
anglefunc = interp1([2,7,17,25],angleraw,1:20,'spline');

%anglerange = anglefunc(pla);

if pla <= 20
    anglerange = anglefunc(pla);
else
    anglerange = 178;
end

switch pla
    case 2
        anglerange = 90;
    case 3
        anglerange = 135;
    case 4
        anglerange = 155;
end

if nargin >= 5
    anglerange = spreadangle;
end

if anglerange == 360
    anglerange = 360-(360/(pla+1));
end


dOm = anglerange*deg/(pla-1);
ddu = 360*deg/spp;
du = 2/3*ddu; %du = 0;

NSAT = pla*spp;

pri = PRIMARY_STATE;

a0 = pri(1);
Om0 = pri(4);
in0 = pri(5);
u0 = pri(6);

STATE510 = zeros(NSAT,6);

STATE510(:,1) = a0;
STATE510(:,2) = 0;
STATE510(:,3) = nan;
STATE510(:,5) = in0;

sats = int16(1:NSAT);

%% ASSIGNING RAAN

STATE510(:,4) = Om0 + double(idivide(sats-1,spp))*dOm;

%% ASSIGNING NEIGHBOUR SHIFT

for i = 1:pla
    
    if mod(i,3) == 1
        STATE510((spp*(i-1)+1):spp*i,6) = u0;
    elseif mod(i,3) == 2
        STATE510((spp*(i-1)+1):spp*i,6) = u0 + du;
    else 
        STATE510((spp*(i-1)+1):spp*i,6) = u0 + du/2;
    end
    
%     switch i
%         case {1,4,7,10,13,16,19,22,25,28,31}
%             STATE510((spp*(i-1)+1):spp*i,6) = u0;
%         case {2,5,8,11,14,17,20,23,26,29,32}
%             STATE510((spp*(i-1)+1):spp*i,6) = u0 + du;
%         case {3,6,9,12,15,18,21,24,27,30,33}
%             STATE510((spp*(i-1)+1):spp*i,6) = u0 + du/2;
%     end
end

%% ASSIGNING PHASE

STATE510(:,6) = STATE510(:,6) + double(mod(sats-1,spp))'*ddu;
STATE510(:,6) = mod(STATE510(:,6),2*pi);

%% DEFINING CARTESIAN

[XYZ,VVV] = math.randv(STATE510(:,1),STATE510(:,2),STATE510(:,5),STATE510(:,4),zeros(NSAT,1),STATE510(:,6));
XYZ = XYZ'; VVV = VVV';

CART510 = [XYZ,VVV];

%% DEFINING GEODETIC

LAT = asin(XYZ(:,3)./sqrt(XYZ(:,1).^2+XYZ(:,2).^2+XYZ(:,3).^2));
THETA = atan2(XYZ(:,2),XYZ(:,1));

THETAG = astro.gstime(jd);

LON = THETA - THETAG;
LON = mod(LON,2*pi);
LON(LON > pi) = LON(LON > pi) - 2*pi;

LATLON510 = [LAT,LON]/deg;

end

