function [ output_args ] = gettrafficsat( a,e,om,Om,in,u, R, jd, deltat, map, mapdata, BIGASSMASK, netdata )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

dOm = netdata{1};
du = netdata{2};
ddu = netdata{3};

[ sectionmask, plotcover, lat, lon ] = getsect( a,e,om,Om,in,u, R, mapdata, jd, dOm, du, ddu );

tz = -timezone(lon,'degrees');

utchr = mod(jd-0.5,1)*24;
localhr = utchr + tz;
localhr = mod(localhr,24);

[ vdata, ddata, sdata, SESTABLE ] = gettraffic( lat, lon, R, localhr, localhr+deltat/3600, map, mapdata, BIGASSMASK, sectionmask, 1 );

for i = 2:7
    FIS = plotcover(:,i,1);
    LAS = plotcover(:,i,2);
    plot(LAS,FIS,'b', 'LineWidth',2, 'LineStyle','--');
end

end

