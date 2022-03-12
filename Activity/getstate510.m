function [ STATE510, CART510, LATLON510 ] = getstate510( PRIMARY_STATE, jd, netdata )
%GET ORBITAL POSITIONS OF 510 SATS
%   PLANE #9 is PRIMARY
%   RETURNS CIRCULAT POSITIONS, IGNORES ECCENTRICITY

deg = pi/180;

dOm = netdata{1};
du = netdata{2};
ddu = netdata{3};

pri = PRIMARY_STATE;

a0 = pri(1);
Om0 = pri(4);
in0 = pri(5);
u0 = pri(6);

STATE510 = zeros(510,6);

STATE510(:,1) = a0;
STATE510(:,2) = 0; 
STATE510(:,3) = nan;
STATE510(:,5) = in0;

sats = int16(1:510);

%% ASSIGNING RAAN

STATE510(:,4) = Om0 - 8*dOm + double(idivide(sats-1,30))*dOm;

%% ASSIGNING NEIGHBOUR SHIFT

for i = 1:17
    switch i
        case {1,4,7,10,13,16}
            STATE510((30*(i-1)+1):30*i,6) = u0 + du;
        case {2,5,8,11,14,17}
            STATE510((30*(i-1)+1):30*i,6) = u0 + du/2;
        case {3,6,9,12,15}
            STATE510((30*(i-1)+1):30*i,6) = u0;
    end
end

%% ASSIGNING PHASE

STATE510(:,6) = STATE510(:,6) + double(mod(sats-1,30))'*ddu;
STATE510(:,6) = mod(STATE510(:,6),2*pi);

%% DEFINING CARTESIAN

[XYZ,VVV] = randv(STATE510(:,1),STATE510(:,2),STATE510(:,5),STATE510(:,4),zeros(510,1),STATE510(:,6));
XYZ = XYZ'; VVV = VVV';

CART510 = [XYZ,VVV];

%% DEFINING GEODETIC

LAT = asin(XYZ(:,3)./sqrt(XYZ(:,1).^2+XYZ(:,2).^2+XYZ(:,3).^2));
THETA = atan2(XYZ(:,2),XYZ(:,1));

THETAG = gstime(jd);

LON = THETA - THETAG;
LON = mod(LON,2*pi);
LON(LON > pi) = LON(LON > pi) - 2*pi;

LATLON510 = [LAT,LON]/deg;

end

