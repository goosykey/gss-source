function [ vdata, ddata, sdata, SESTABLE ] = gettrafficsubs( t, SUBSAT, map, mapdata, BIGASSMASK, netdata, utcmodif )
%GETTRAFFIC by click on SUBSCRIBERS PLOT FROM SIMULATION RESULT
%   TICK IS 1 SECOND
%   specify NEGATIVE HR for UTC TIME

if nargin < 7
    utcmodif = 0;
end

[HR,~] = ginput(1);

[~, mintc] = min(abs(t-HR*3600));

lat = SUBSAT(mintc,1);
lon = SUBSAT(mintc,2);

HR = mod(HR,24) + utcmodif;

[ vdata, ddata, sdata, SESTABLE ] = gettraffic( [lat lon 750], -HR, map, mapdata, BIGASSMASK, netdata, 1 );


end

