function [ isinside ] = belongstorect( FI, LA, fi1, fi2, la1, la2 )
%IF POINTS BELONG TO GEODETIC RECTANGLE (Returns boolean)
%   Detailed explanation goes here

if fi1 < fi2
    FI(FI(:)>fi2) = 0;
    FI(FI(:)<fi1) = 0;
else
    FI(FI(:)>fi1) = 0;
    FI(FI(:)<fi2) = 0;
end

if la1 < la2
    LA(LA(:)>la2) = 0;
    LA(LA(:)<la1) = 0;
else
    LA(LA(:)>la1) = 0;
    LA(LA(:)<la2) = 0;
end

isinside = boolean (FI.*LA);


end

