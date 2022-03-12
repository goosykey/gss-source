function [ ifsees, distances, alphas ] = satsees( Recf, latlons, elevcrit )
%IF SATELLITE SEES POINT(ARRAY)
%   Use to determine whether sat sees an array of points
%   INPUT:
%       Recf - cartesian position of the satellite in ECF, kilometres
%       latlons - (Nx2) array of point coordinates, deg
%       elevcrit - critical elevation in degrees (optionaldefault 0)
% 
%   OUTPUT: 
%       ifsees - (Nx1) boolean - if sees 
%       distances - km
%       alphas - degrees
% 

if nargin < 3
    elevcrit = 0;
end

N = length(latlons(:,1));

R0 = 6378.137;

R = sqrt(sum(Recf.^2));

alt = norm(Recf) - R0;
alphacrit = asind(sind(90 + elevcrit) * R0 / (R0 + alt));

POINTS = zeros(N,3);
POINTS(:,3) = R0.*sind(latlons(:,1));
POINTS(:,2) = R0.*cosd(latlons(:,1)).*sind(latlons(:,2));
POINTS(:,1) = R0.*cosd(latlons(:,1)).*cosd(latlons(:,2));

RecfN = repmat(Recf,N,1);

R01s = POINTS - RecfN;
distances = sqrt(sum(R01s.*R01s,2));

cosalphas = dot(-RecfN,R01s,2) ./ R ./ sqrt(sum(R01s.^2,2));
alphas = acosd(cosalphas);

ifsees = (alphas < alphacrit & distances < 3000);



end

