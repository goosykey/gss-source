function [ Ndes,Nreal, EPHrnd ] = accesshere( lat, lon, jd, netdata, map, mapdata, BIGASSMASK )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

R0 = 6378.14; deg = pi/180; mu = 398600.44;

a = R0+600; e = 0; om = 0; in = 97.79*deg;

latlim = mapdata{1}; lonlim = mapdata{2};
latstep = mapdata{3}; lonstep = mapdata{4};
M = mapdata{5}; N = mapdata{6}; cliffname = mapdata{7};

lats = 90:-latstep:-90;
lons = -180:lonstep:180;
[latmap, lonmap] = ndgrid(lats,lons);

di = R0*deg*distance(latmap,lonmap,lat, lon);

diboo = di < 750;

t = 0:60:86400; t=t'; NN = numel(t);
EPHEMERIS = [t,zeros(NN,6)];
EPHEMERIS(:,2) = a;
EPHEMERIS(:,3) = e;
EPHEMERIS(:,4) = om;
EPHEMERIS(:,6) = in;
fprintf('Generating ephemeris...\n');

for q = 0:60:86400
    fprintf(' > %g', q);
    [lonsub,latsub] = pinky(lons,lats,diboo,3);
    
    Recf = [R0*cosd(latsub)*cosd(lonsub) R0*cosd(latsub)*sind(lonsub)...
        R0*sind(latsub)]'*a/R0;
    
    Reci = R3(-gstime(jd+q/86400))*Recf;
    v = sqrt(mu/a); x = Reci(1); y = Reci(2); z = Reci(3); r = a;
    
    syms vx vy vz
    [sol_vx, sol_vy, sol_vz] = vpasolve([vx^2 + vy^2 + vz^2 == v^2,...
        x*vx+y*vy+z*vz==0, x*vy-y*vx==r*v*cos(in)], [vx, vy, vz]);
    
    
    [ ~,~,~,Om,~,u ] = xyz2kepler( x,y,z,double(sol_vx(1)), double(sol_vy(1)), double(sol_vz(1)) );
    
    EPHEMERIS(t==q,5) = Om;
    EPHEMERIS(t==q,7) = u;
    
    
end
fprintf('\n');

[ Nreal, Ndes, ~ ] = CYCLO_SIM( EPHEMERIS, netdata, jd, map, mapdata, BIGASSMASK  ); % get P_out from session distributions

EPHrnd = EPHEMERIS;

end

