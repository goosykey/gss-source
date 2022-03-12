function [ map, mapdata ] = readmap( filename )
%GET MAP AND MAPDATA FROM .asc population file
%   Detailed explanation goes here

[~ , val] = textread(filename,'%s %f', 6, 'delimiter',' ','headerlines',0,'emptyvalue',0);
latlim = [0,0];
lonlim = [0,0];

ncols = val(1);
nrows = val(2);
cellsize = val(5);

lonlim(1) = val(3);
latlim(2) = val(4);

lonlim(2) = round(lonlim(1) + ncols*cellsize);
latlim(1) = round(latlim(2) + nrows*cellsize);

latstep = cellsize;
lonstep = cellsize;

fprintf('Reading map file... \n');

map = dlmread(filename,' ',6,0); 
map(:,length(map(1,:))) = [];

[M,N] = size(map);

mapdata = {latlim, lonlim, latstep, lonstep, M, N, 'cliffs.csv'};


end

