function [ BIGASSMASK ] = intmaskbatch( map, mapdata, cdb, uidb )
%USAGE INTENSITY MASK INITIALIZATION
%   Detailed explanation goes here

%% INITIALIZE MAP PARAMETERS

latlim = mapdata{1}; lonlim = mapdata{2};
latstep = mapdata{3}; lonstep = mapdata{4};

if latlim(1)>latlim(2)
    latstep = -latstep;
end

LATS = latlim(1):latstep:latlim(2);
LONS = lonlim(1):lonstep:lonlim(2);

[latmap, lonmap] = ndgrid(LATS,LONS);

latL = length(LATS);
lonL = length(LONS);

Vmask = zeros(latL,lonL);
Dmask = zeros(latL,lonL);
Smask = zeros(latL,lonL);
Rmask = zeros(latL,lonL);
Pmask = zeros(latL,lonL);
SPmask= zeros(latL,lonL);
Nmask = zeros(latL,lonL);
Zmask = zeros(latL,lonL);

latlonrr = cell2mat(cdb(:,4:7));
ncities = length(latlonrr(:,1));

latlonrr = [latlonrr,(1:ncities)'];

latlonrr = sortrows(latlonrr,1);

%% FORM OCEAN MASK

tic;
fprintf('Forming ocean mask...\n');
%Omask = map*0;
Omask = oceanALL(latmap,lonmap,40);
T = toc;
fprintf('Mask formed in %3.2f sec.\n',T);


%% MAIN CYCLE - REVERSE GEOCODE EACH POINT

% figure('Name','LATITUDE SCANNING','NumberTitle','on'); hold on;
% area(LATS, LATS./LATS-1);

parpool(2);
fprintf('Latitude scanning begins... \n');

hon = 450; % honesty of 350km

urbanc = 0.04;
ruralc = 0.25;
dstntc = 0.8;
%global temp;

parfor q = 1:latL
    latcur = LATS(q);
    fprintf('Scanning %02.3f deg... \n',latcur);
    latcropped = latlonrr(abs(latlonrr(:,1)-latcur)<hon/121,:);
    latcropped = sortrows(latcropped,2);
    %disp(numel(latcropped(:,1)));
    
    [range, index, zone, pencoef] = findcitybatch( latcur, LONS, latcropped, hon );
    
    isocode = cell(numel(LONS),1);
    isocode(~isnan(index)) = cdb(index(~isnan(index)),1);
    isocode(isnan(index)) = {'AQ'};
    %temp = isocode;
    
    denshere = map(q,:)/cellarea(q,mapdata); % POP DENSITY IN THIS ROW
    Rmask(q,:) = range;
    
    pencoef(denshere > 750 & zone < 2) = 0.2;
    zone(denshere > 750 & zone < 2) = 2;
    pencoef(denshere > 150 & zone < 1) = 1;
    zone(denshere > 150 & zone < 1) = 1;
    
    zone(denshere < 10 & Omask(q,:)) = zone(denshere < 10 & Omask(q,:)) + 4;
    
    [ VUI, DUI, SUI, PEN, SPEN, CountryNo ] = getuisbatch( isocode, uidb );
    
    Pmask(q,:) = PEN'.*pencoef.*(1+(zone>2)*2);
    SPmask(q,:) = SPEN'.*pencoef.*(1+(zone>2)*2);
    Nmask(q,:) = CountryNo';
    
    Zmask(q,:) = zone;
    
    zonerow = (zone==0)*dstntc + (zone==1)*ruralc + (zone==2)*urbanc...
        + (zone==4)*2 + (zone==5)*1.5 + (zone==6)*1.25;
    
    Vmask(q,:) = VUI'.*zonerow * 1.0;
    Dmask(q,:) = DUI'.*zonerow * 1.2;
    Smask(q,:) = SUI'.*zonerow * 2.0;
    
end % endfor lat

delete(gcp('nocreate'));

Vmask = int16(Vmask*100);
Dmask = int16(Dmask*100);
Smask = int16(Smask*100);
Nmask = int16(Nmask);
Zmask = int8(Zmask);

%% ADJUST PENETRATION MASK (NORMALIZE)

[ PmaskAd, SPmaskAd ] = penadjust( Pmask, SPmask, Nmask, map, uidb );
%why do this inside the function - because I can!

%% PUT EVERYTHING IN A CELL ARRAY BECAUSE WHY THE FUCK NOT

BIGASSMASK = {Vmask, Dmask, Smask, Rmask, Omask, PmaskAd, SPmaskAd, Nmask, Zmask};


end

