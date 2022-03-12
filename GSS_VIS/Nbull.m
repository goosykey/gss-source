function [ Nmean, Nmax, Nmin ] = Nbull( CARTN, num )
%NBULL Cross-coverage analysis for SPACE Y
%   Detailed explanation goes here

if nargin < 2
    num = 10000;
end

R0 = 6378.14;
anglecrit = 50;


[Nsat, ~] = size(CARTN);
Ncol = zeros(num,1);

for i = 1:num
    lat = 2*53*rand - 53;
    lon = 360 * rand;
    X = R0 * cosd(lat) * cosd(lon);
    Y = R0 * cosd(lat) * sind(lon);
    Z = R0 * sind(lat);
    XYZ = repmat([X Y Z],[Nsat,1]);
    R1s = CARTN-XYZ;
    DOT = dot(XYZ,R1s,2);
    XYZabs = sqrt(sum(XYZ.^2,2));
    R1sabs = sqrt(sum(R1s.^2,2));
    angles = acosd(DOT./XYZabs./R1sabs);
    Ncol(i) = sum(angles < anglecrit);
end

Nmean = mean(Ncol);
Nmax = max(Ncol);
Nmin = min(Ncol);

end

