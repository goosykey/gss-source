function [ Lmean ] = meanvis595( latdeg )
%MEANVIS168 Summary of this function goes here
%   Detailed explanation goes here

if any(latdeg < -90 | latdeg > 90)
    error ('gtfo');
end

Lmeanraw = [1.086 1.112 1.241 1.548 2.246 5.686 4.3];

Lmeanfunc = interp1(0:15:90,Lmeanraw,0:1:90,'spline');

Lmean = Lmeanfunc(round(abs(latdeg))+1);

Lmean = Lmean';

end

