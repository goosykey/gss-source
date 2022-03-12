function [cdbmodified, poplatlon] = popmodify(census, dspot)
%POPMODIFY Modify populations in spot
%   Detailed explanation goes here

fprintf('Modifying population in spot...\n');

R0 = 6378.14;
deg = pi/180;

%% INITIALIZE MAP PARAMETERS

[ map, mapdata ] = readmap( 'glp15ag.asc' );
cdb = readcities('cities1000.txt', census);

cdb(strcmp(cdb(:,2) , {'City of London'}),:)=[]; %% REMOVE THIS SHIT

poplatlon = cell2mat(cdb(:,3:5));
Nc = length(poplatlon(:,1));
poplatlon = [poplatlon,(1:Nc)'];

latlim = mapdata{1}; lonlim = mapdata{2};
latstep = mapdata{3}; lonstep = mapdata{4};

if latlim(1)>latlim(2)
    latstep = -latstep;
end

LATS = latlim(1):latstep:latlim(2);
LONS = lonlim(1):lonstep:lonlim(2);

[latmap, lonmap] = ndgrid(LATS,LONS);

%% MODIFY POPULATION - POP IN THE SPOT

for q = 1:Nc
    LAT = poplatlon(q,2);
    LON = poplatlon(q,3);
    DISTMAT = R0 .* deg .* distance(LAT,LON,latmap,lonmap);
    maphere = map;
    maphere(DISTMAT > dspot/2) = 0;
    pophere = sum(sum(maphere));
    poplatlon(q,1) = pophere;
end

fprintf('done\n');

%% MODIFY CDB 

cdbmodified = cdb;

for q = 1:Nc
    cdbmodified{q,3} = poplatlon(q,1);
end

end

