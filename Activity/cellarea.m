function [ S ] = cellarea( i,mapdata )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

latlim = mapdata{1}; latstep = mapdata{3};

lat1 = latlim(1);

lat = lat1 - latstep * (i-1);

st1 = 6378.14*pi/180*abs(latstep);

S = st1^2*cosd(lat);


end

