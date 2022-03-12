function [ plotcover ] = getcoverage( lat, lon, h, R, netdata, crop )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

R0 =6378.14;
deg = pi/180;
mu = 398600.44;

dOm = netdata{1};
du = netdata{2};
ddu = netdata{3};

a = R0 + h;
e = 0; om = 0;
in = 97.792*deg;

Recf = [R0*cosd(lat)*cosd(lon) R0*cosd(lat)*sind(lon)...
    R0*sind(lat)]'*a/R0;

Reci = Recf; % don't care much;
v = sqrt(mu/a); x = Reci(1); y = Reci(2); z = Reci(3); r = a;

syms vx vy vz
[sol_vx, sol_vy, sol_vz] = vpasolve([vx^2 + vy^2 + vz^2 == v^2,...
    x*vx+y*vy+z*vz==0, x*vy-y*vx==r*v*cos(in)], [vx, vy, vz]);
[ ~,~,~,Om,~,u ] = xyz2kepler( x,y,z,double(sol_vx(1)), double(sol_vy(1)), double(sol_vz(1)) );

% WARNING: 'Om' IN THIS FUNCTION IS ACTUALLY LAN (LON.ASC.NODE) !!!

A = [a a a a a a a];
E = [e e e e e e e];
PERIG = [om om om om om om om];
OMEGA = [Om Om Om+dOm Om+dOm Om Om-dOm Om-dOm];
IN = [in in in in in in in];
U = [u u+ddu u+du u+du-ddu u-ddu u-du u-du+ddu];

[Rxyz,~] = randv(A,E,IN,OMEGA,PERIG,U-PERIG); % getting in ECF

t = 0:1:360;
plotcover = zeros(length(t),7,2);

lat1 = lat; lon1 = lon;

for i = 1:7
    Ri = Rxyz(:,i); % already in ECF
    ri = norm(Ri);
    lat = asin(Ri(3)/ri)/deg;
    lon = atan2(Ri(2),Ri(1))/deg;
    gamma = atan(R/R0);
    zg = R0*cos(gamma);
    xg = R0*sin(gamma)*cosd(t);
    yg = R0*sin(gamma)*sind(t);
    GG = zeros(3,length(t));
    GG(1,:) = xg; GG(2,:) = yg; GG(3,:) = zg;
    XYZ = R3(pi/2-lon*pi/180)*R1(pi/2-lat*pi/180)*GG;
    FIS = asin(XYZ(3,:)/R0)*180/pi;
    LAS = atan2(XYZ(2,:),XYZ(1,:))*180/pi;
    LAS(abs(LAS(:)-lon)>180) = LAS(abs(LAS(:)-lon)>180) + 360;
    LAS((LAS(:)-lon)>180) = LAS((LAS(:)-lon)>180) - 720; 
    if crop
        DIST = R0*deg*distance(FIS, LAS, lat1, lon1);
        LAS(DIST > R) = nan;
        FIS(DIST > R) = nan;
    end
    plotcover(:,i,1) = FIS;
    plotcover(:,i,2) = LAS;
end

end

