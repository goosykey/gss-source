function [ visible, alphas, elevations, distances ] = ifvisible ( point, satarray, elevcrit )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if nargin < 3
    elevcrit = 0;
end

if numel(point) ~= 3
    error('gtfo point');
end

point = [point(1) point(2) point(3)];

[M,N] = size(satarray);

if not(any([M N] == 3))
    error('gtfo sats');
end

if N ~= 3 && M == 3
    satarray = satarray';
end

[M,~] = size(satarray);

R0 = 6378.137;
alt = norm(satarray(1,:)) - R0;
alphacrit = asind(sind(90 + elevcrit) * R0 / (R0 + alt));

R0ECIN = repmat(point,M,1);
R01N = satarray - R0ECIN;
dists = sqrt(sum(R01N.*R01N,2));
rads = sqrt(sum(satarray.*satarray,2));
DOTN = dot(R01N,satarray,2);
cosalphas = DOTN ./ dists ./rads;
alphas = acosd(cosalphas);
visible = ( alphas < alphacrit & dot(R01N,R0ECIN,2) > 0 );

n = point ./ R0;
nn = repmat(n,M,1);
elevations = 90 - acosd(dot(nn,R01N,2)./sqrt(sum(R01N.^2,2)));
distances = dists;


end

