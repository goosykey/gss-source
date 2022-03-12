function [ BIGASSMASK ] = intmasknew( map, mapdata, cdb, uidb )
%USAGE INTENSITY MASK INITIALIZATION
%   20-05-2016
%   Alexander Kharlan   |   Yaliny
%   OUTPUT:
%   BIGASSMASK = {Vmask, Dmask, Smask, Rmask, Omask, PmaskAd, SPmaskAd, Nmask, Zmask};
%   Vmask, Dmask, Smask - usage intensity
%   Rmask - map of ranges to nearest cities
%   Omask - boolean mask defining whether the cell is ocean or not
%   PmaskAd, SPmaskAd - subscriber penetration mask (percentage of
%       subscribers in each cell)
%   Nmask - mask of country IDs based on provided uidb table
%   Zmask - zone type mask (2 - urban, 1 - rural, 0 - distant, 4..6 -
%       ocean)
%   
%   INPUTS:
%   map, mapdata - formed by 'readmap.m'
%   cdb - formed by 'readcities.m'
%   uidb - formed by 'readui.m'
%   
%   WARNING: this version of INTMASK uses parallel calculations thus
%       possibly causing severe lags and RAM overflow. To avoid this, use
%       'intmaskone.m' instead

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
Omask = oceanALL(latmap,lonmap,40);
T = toc;
fprintf('Mask formed in %3.2f sec.\n',T);


%% MAIN CYCLE - REVERSE GEOCODE EACH POINT

% figure('Name','LATITUDE SCANNING','NumberTitle','on'); hold on;
% area(LATS, LATS./LATS-1);

parpool(3);
fprintf('Parallel latitude scanning begins... \n');

hon = 450; % honesty of 350km
urbanc = 0.04;
ruralc = 0.20;
dstntc = 0.8;

parfor q = 1:latL
    latcur = LATS(q);
    fprintf('Scanning %g deg... \n',latcur);
    latcropped = latlonrr(abs(latlonrr(:,1)-latcur)<hon/121,:);
    %latcropped = sortrows(latcropped,2);
    for w = 1:lonL
        loncur = LONS(w);
        
        loncropped = latcropped(abs(latcropped(:,2)-loncur)<hon/121/cosd(latcur),:);
        
        [range, index, zone, pencoef] = findcity( latcur, loncur, loncropped, hon );
        if isnan(index)
            isocode = 'AQ';
        else
            isocode = cdb{index,1};
        end
        denshere = map(q,w)/cellarea(q,mapdata); % POP DENSITY IN THIS CELL
        Rmask(q,w) = range;
        if denshere > 750 && zone < 2  % CONVERT TO URBAN
            pencoef = 0.2;
            zone = 2;
        elseif denshere > 150 && zone < 1 % CONVERT TO RURAL
            pencoef = 1;
            zone = 1;
        end
        if denshere < 10 && Omask(q,w)
            zone = zone + 4;
        end
        [ VUI, DUI, SUI, PEN, SPEN, CountryNo ] = getuis( isocode, uidb );
        Pmask(q,w) = PEN*pencoef*(1+(zone>2)*2);
        SPmask(q,w) = SPEN*pencoef*(1+(zone>2)*2);
        Nmask(q,w) = CountryNo;       
        
        Zmask(q,w) = zone;
        % SWITCH ZONE BEGINS
        switch zone
            case 0 % DISTANT
                Vmask(q,w) = VUI*dstntc * 1.0;
                Dmask(q,w) = DUI*dstntc * 1.0;
                Smask(q,w) = SUI*dstntc * 2.5;
            case 1 % RURAL
                Vmask(q,w) = ruralc*VUI;
                Dmask(q,w) = ruralc*DUI * 0.3;
                Smask(q,w) = ruralc*SUI * 2.0;
            case 2 % URBAN
                Vmask(q,w) = urbanc*VUI;
                Dmask(q,w) = urbanc*DUI * 0.3;
                Smask(q,w) = urbanc*SUI;
            case 4 % DISTANT OCEAN
                Vmask(q,w) = 2*VUI;
                Dmask(q,w) = 3*DUI;
                Smask(q,w) = 4*SUI;
            case 5 % "RURAL" OCEAN
                Vmask(q,w) = 1.5*VUI;
                Dmask(q,w) = 1.5*DUI;
                Smask(q,w) = 2*SUI;
            case 6 % "URBAN" OCEAN
                Vmask(q,w) = 1.2*VUI;
                Dmask(q,w) = 0.7*DUI;
                Smask(q,w) = 2*SUI;
            otherwise % ???
                Vmask(q,w) = 2;
                Dmask(q,w) = 2;
                Smask(q,w) = 5;
        end % end switch zone
        % SWITCH ZONE ENDS
    end % endfor long
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

