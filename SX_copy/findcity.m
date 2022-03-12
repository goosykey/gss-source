function [ range, index, zone, pencoef ] = findcity( lat, lon, latlonrr, honesty )
%REVERSE GEOCODE MAIN SUBFUNCTION
%     INPUT:
%         lat 	 : point latitude
%         lon 	 : point longitude
%         latlonrr : city database lat/lon dump + two zone radiuses for each (Nx5 vector)
%                       (see "COMPUTE ZONAL MAP" in 'readcities.m')
%                       latlonrr is cropped in 'intmasknew.m' to speed up
%                       the algotithm
%         honesty : <optional> ranges (km) bigger than this will be considered Inf.
%     OUTPUT:
%         range	 : distance (some kind of) to the nearest settlement
%         index	 : index of the settlement in the provided 'latlonrr' col.5
%         zone	 : type of the zone:
%                 0 : distant, 4 : distant, ocean
%                 1 : rural,   5 : rural, ocean (shore)
%                 2 : urban    6 : urban, ocean (city coast)
%         pencoef: subs penetr in this point to country average (ratio)

%% ARGUMENT OPERATIONS

if nargin < 4
    honesty = 300;
end

if isempty(latlonrr)
    range = inf;
    index = NaN;
    zone = 0;
    pencoef = 10;
    return
end

%% INITIALIZE MAP PARAMETERS

latsr = latlonrr(:,1)*pi/180; latr = lat*pi/180;
lonsr = latlonrr(:,2)*pi/180; lonr = lon*pi/180;


%% TRY TO RETURN ON BIG RANGES

% PSEUDO-RANGES - MERCATOR RECTANGULAR DISTANCE
% pranges = 6378.14*pi/180 * sqrt((lats-lat).^2 + ((lons-lon).*cosd(lat)).^2);
% 
% if ~any(pranges < honesty*1.2)
%     range = inf;
%     index = NaN;
%     zone = 0;
%     pencoef = 4;
%     return
% end

% Great Circle Distances - TRUE RANGES
% GCD option finds arc distance, takes MUCH longer

%ranges = 6378.14*pi/180*distance(latlonrr(:,1),latlonrr(:,2),lat,lon);

% SEMI-TRUE RANGES - CARTESIAN DISTANCE THROUGH CORE

ranges = 6378.14*sqrt((cos(latsr).*cos(lonsr) - cos(latr).*cos(lonr)).^2 ...
    + (cos(latsr).*sin(lonsr) - cos(latr).*sin(lonr)).^2 ...
    + (sin(latsr) - sin(latr)).^2);

if ~any(ranges < honesty)
    range = inf;
    index = NaN;
    zone = 0;
    pencoef = 4;
    return
end

%% DETERMINE ZONE CATEGORIES

% heaviside takes a little longer
%zonecat = heaviside(latlonrr(:,3)-ranges) + heaviside(latlonrr(:,4)-ranges);

%tic
urban = latlonrr(:,3)-ranges;
urban(urban(:)<=0) = 0;
urban(urban(:)>0) = 1;

rural = latlonrr(:,4)-ranges;
rural(rural(:)<=0) = 0;
rural(rural(:)>0) = 1;

zonecat = urban + rural;


if any(zonecat==2)
    indexes = (1:length(ranges))';
    bigarray = [indexes, ranges, zonecat];
    deux = bigarray(bigarray(:,3)==2,:);
    [range,mini] = min(deux(:,2));
    index = deux(mini,1);
    pencoef = 0.2;
elseif any(zonecat==1)
    indexes = (1:length(ranges))';
    bigarray = [indexes, ranges, zonecat];
    un = bigarray(bigarray(:,3)==1,:);
    [range,mini] = min(un(:,2));
    index = un(mini,1);
    R1 = latlonrr(index,3);
    R2 = latlonrr(index,4);
    pencoef = 0.2 + (range - R1)*1.3/(R2 - R1);
else
    [range,index] = min(ranges);
    R1 = latlonrr(index,4);
    R2 = latlonrr(index,4)*3;
    pencoef = 1.5 + (range - R1)*2.5/(R2 - R1);
    if pencoef > 4
        pencoef = 4;
    end
end

zone = zonecat(index);

index = latlonrr(index,5);


end

