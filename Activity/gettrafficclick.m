function [ vdata, ddata, sdata, SESTABLE ] = gettrafficclick( HR, map, mapdata, BIGASSMASK, netdata )
%GETTRAFFIC by google map click
%   TICK IS 1 SECOND
%   specify NEGATIVE HR for UTC TIME

if  isempty(get(gcf,'CurrentAxes'))
    plot_google_map;
end
[lon,lat] = ginput(1);

[ vdata, ddata, sdata, SESTABLE ] = gettraffic( [lat lon 1060], HR, map, mapdata, BIGASSMASK, netdata, 33 );


end

