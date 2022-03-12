function [ range, index, zone, pencoef ] = findcitybatch( lat, lonrow, latlonrr, honesty )
%REVERSE GEOCODE MAIN SUBFUNCTION
%     INPUT:
%         lat 	 : point latitude
%         lon 	 : point longitude
%         latlonrr : city database lat/lon dump + two zone radiuses for each (Nx4 vector)
%                       (see "COMPUTE ZONAL MAP" in 'readcities.m')
%         honesty : <optional> ranges (km) bigger than this will be considered Inf.
%     OUTPUT:
%         range	 : distance (some kind of) to the nearest settlement
%         index	 : index of the settlement in the provided 'latlonrr'
%         zone	 : type of the zone:
%                 0 : distant, 4 : distant, ocean
%                 1 : rural,   5 : rural, ocean (shore)
%                 2 : urban    6 : urban, ocean (city coast)
%         pencoef: subs penetr in this point to country average (ratio)

%% ARGUMENT OPERATIONS

if isempty(latlonrr)
    range = ones(1,numel(lonrow))*inf;
    index = zeros(1,numel(lonrow))*nan;
    zone = zeros(1,numel(lonrow));
    pencoef = ones(1,numel(lonrow))*4;
    return
end

%% PREALLOCATE OUTPUT

range = zeros(1,numel(lonrow));
index = zeros(1,numel(lonrow));
pencoef = zeros(1,numel(lonrow));

%% MATRIX MANIPULATIONS & RESHAPES

nc = numel(latlonrr(:,1));
nd = numel(lonrow);

latcity = repmat(latlonrr(:,1),[nd,1]);
loncity = repmat(latlonrr(:,2),[nd,1]);

latdot = lat*ones(1,nd);
latdot = repmat(latdot,[nc,1]);
latdot = reshape(latdot,[nc*nd,1]);

londot = repmat(lonrow,[nc,1]);
londot = reshape(londot,[nc*nd,1]);

% LLLL = [latcity, loncity, latdot,londot];
% LLLL(abs(loncity-londot) < honesty / 121 / cosd (lat),:) = nan;
% 
% latcity = LLLL(:,1); loncity = LLLL(:,2);
% latdot = LLLL(:,3); londot = LLLL(:,4);
% clear LLLL;

RANGEMAT = single(6378.14*sqrt((cosd(latcity).*cosd(loncity) - cosd(latdot).*cosd(londot)).^2 ...
    + (cosd(latcity).*sind(loncity) - cosd(latdot).*sind(londot)).^2 ...
    + (sind(latcity) - sind(latdot)).^2));

RANGEMAT = reshape(RANGEMAT,[nc,nd]);

clear latcity loncity latdot londot

RuMAT = repmat(latlonrr(:,3),[1,nd]);
RrMAT = repmat(latlonrr(:,4),[1,nd]);



%% DETERMINE ZONE CATEGORIES

urban = RuMAT-RANGEMAT;
urban(urban(:)<=0) = 0;
urban(urban(:)>0) = 1;

rural = RrMAT-RANGEMAT;
rural(rural(:)<=0) = 0;
rural(rural(:)>0) = 1;

zonecat = urban + rural;

isurban = sum(zonecat==2,1) > 0;
isrural = sum(zonecat==1,1) > 0 - isurban;

zone = int8(isurban/2 + isrural);

clear RrMAT RuMAT



rangesur = RANGEMAT.*(zonecat==2); rangesur(rangesur==0) = nan;
[minr,mini] = min(rangesur,[],1);
range(zone==2) = minr(zone==2);
index(zone==2) = mini(zone==2);
clear rangesur
%fprintf('%g - ',max(sum(zonecat==2)));

rangesru = RANGEMAT.*(zonecat==1); rangesru(rangesru==0) = nan;
[minr,mini] = min(rangesru,[],1);
range(zone==1) = minr(zone==1);
index(zone==1) = mini(zone==1);
clear rangesru
%fprintf('%g \n ',max(sum(zonecat==1)));

rangesdi = RANGEMAT.*(zonecat==0); rangesdi(rangesdi==0) = nan;
[minr,mini] = min(rangesdi,[],1);
range(zone==0) = minr(zone==0);
index(zone==0) = mini(zone==0);
clear rangesdi

pencoef(zone==2) = 0.2;



R1 = latlonrr(index,3)';
R2 = latlonrr(index,4)';
pencoef(zone==1) = 0.2 + (range(zone==1) - R1(zone==1))*1.3/(R2(zone==1) - R1(zone==1));

R1 = latlonrr(index,4)';
R2 = latlonrr(index,4)'*3;
pencoef(zone==0) = 0.2 + (range(zone==0) - R1(zone==0))*1.3/(R2(zone==0) - R1(zone==0));

pencoef(pencoef>4) = 4;

index = latlonrr(index,5)';
range(range > honesty) = inf;
index(range > honesty) = nan;
pencoef(range > honesty) = 4;

end

