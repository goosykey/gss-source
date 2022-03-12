function [  ] = lookupclick( mapdata, cdb )
%HELP FUNCTION TO ANALYZE MASKS
%   Detailed explanation goes here
%     INPUT
%         lat, lon    : coordinates of the desired point on a map
%         mapdata : data for the given map
%         cdb     : city database (generate thru 'readcities.m')
%     OUTPUT:
%         none

if  isempty(get(gcf,'CurrentAxes'))
    plot_google_map;
end

[lon,lat] = ginput(1);
plot(lon,lat,'dr','MarkerFaceColor','red');

latlim = mapdata{1}; lonlim = mapdata{2};
latstep = mapdata{3}; lonstep = mapdata{4};

lat1 = latlim(1);
lon1 = lonlim(1);

i = (lat1-lat)/latstep + 1; i = round(i);
j = (lon-lon1)/lonstep + 1; j = round(j);


latlonrr = cell2mat(cdb(:,4:7));
ncities = length(latlonrr(:,1));
latlonrr = [latlonrr,(1:ncities)'];

[ range, index, zone ] = findcity( lat, lon, latlonrr, inf );

if isnan(index)
    city = 'NONE';
    country = '';
    popul = NaN;
    latcity = NaN;
    loncity = NaN;
elseif index==0
    fprintf ('The requested point appears to be in the ocean. \n');
    city = 'NONE';
    country = '';
    popul = NaN;
    latcity = NaN;
    loncity = NaN;
else
    city = cdb{index,2};
    country = cdb{index,1};
    popul = cdb{index,3};
    latcity = latlonrr(index,1);
    loncity = latlonrr(index,2);
end




if zone == 2
    zonetype = 'Urban';
elseif zone == 1
    zonetype = 'Rural';
else
    zonetype = 'Distant';
end

fprintf ('Point : %3.3f %3.3f [%g %g on map] \nReference settlement: %s, %s (%3.3f, %3.3f). \n\t Population : %g \n\t Range      : %g km\n\t Zonetype   : %s\n', ...
    lat, lon, i,j, city, country, latcity, loncity, popul, range, zonetype);


end

