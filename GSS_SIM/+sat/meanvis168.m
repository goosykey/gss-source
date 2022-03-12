function [ Lmean ] = meanvis168( latdeg )
%MEANVIS168 Summary of this function goes here
%   Detailed explanation goes here

if any(latdeg < -90 | latdeg > 90)
    error ('gtfo');
end

Lmeanraw = [2.1 2.16 2.43 3.03 4.54 8.58 12.5];

Lmeanfunc = interp1(0:15:90,Lmeanraw,0:1:90,'spline');

Lmean = Lmeanfunc(round(abs(latdeg))+1);

Lmean = Lmean';

end

