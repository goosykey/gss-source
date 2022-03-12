function [ sectionmask, plotcover, lat1, lon1 ] = getsect( a,e,om,Om,in,u, R, mapdata, jd, RAANshift, AOLshift, AOLdelta )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

R0 =6378.14;
deg = pi/180;

A = [a a a a a a a];
E = [e e e e e e e];
PERIG = [om om om om om om om];
OMEGA = [Om Om Om+RAANshift Om+RAANshift Om Om-RAANshift Om-RAANshift];
IN = [in in in in in in in];
U = [u u+AOLdelta u+AOLshift u+AOLshift-AOLdelta u-AOLdelta u-AOLshift u-AOLshift+AOLdelta];

[Rxyz,~] = randv(A,E,IN,OMEGA,PERIG,U-PERIG);

latlim = mapdata{1}; lonlim = mapdata{2};
latstep = mapdata{3}; lonstep = mapdata{4};

lats = latlim(1):-latstep:latlim(2);
lons = lonlim(1):lonstep:lonlim(2);

M = length(lats);
N = length(lons);

[latmap, lonmap] = ndgrid(lats,lons);

sectionmask = zeros(M,N);

t = 0:1:360;
plotcover = zeros(length(t),7,2);

for i = 1:7
    Ri = Rxyz(:,i);
    Ri = R3(gstime(jd))*Ri;
    ri = norm(Ri);
    lat = asin(Ri(3)/ri)/deg;
    lon = atan2(Ri(2),Ri(1))/deg;
    if i==1
        lat1 = lat; lon1 = lon;
    end
    RANGES = 6378.14*pi/180*distance(latmap,lonmap,lat,lon);
    section = zeros(M,N);
    section(RANGES(:)<=R) = 1;
    sectionmask = sectionmask + section;
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
    plotcover(:,i,1) = FIS;
    plotcover(:,i,2) = LAS;
end

sectionmask(sectionmask(:)==0)=inf;
sectionmask = 1./sectionmask;

end

