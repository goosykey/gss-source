function [ rcv ] = rcvcnd( rcv, mapdata, R, addbad, howmany )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if nargin < 2
    error('Not enough input arguments');
elseif nargin < 3
    R = 100;
    addbad = true;
    howmany = 5;
elseif nargin < 4
    addbad = true;
    howmany = 5;
elseif nargin < 5
    howmany = 5;
end

addbad = boolean(addbad);

if  isempty(get(gcf,'CurrentAxes'))
    plot_google_map('maptype' ,'satellite');
end

latlim = mapdata{1}; lonlim = mapdata{2};
latstep = mapdata{3}; lonstep = mapdata{4};

lats = latlim(1):-latstep:latlim(2);
lons = lonlim(1):lonstep:lonlim(2);

[latmap, lonmap] = ndgrid(lats,lons);

[M,N] = size(latmap);

for q = 1:howmany
    [lonc,latc] = ginput(1);
    latm = repmat(latc,[M,N]);
    lonm = repmat(lonc,[M,N]);
    ranges = 6378.14*sqrt((cosd(latmap).*cosd(lonmap) - cosd(latm).*cosd(lonm)).^2 ...
        + (cosd(latmap).*sind(lonmap) - cosd(latm).*sind(lonm)).^2 ...
        + (sind(latmap) - sind(latm)).^2);
    ranges(ranges > R) = nan;
    rcv(~isnan(ranges)) = addbad;
    rgbdot = [addbad, ~addbad, 0];
    plot( lonmap(~isnan(ranges)) , latmap(~isnan(ranges)), 'o', 'Color', rgbdot);    
end


end

