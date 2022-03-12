function [ doy,tod ] = doytod( ST0,t )
%Converts ST0 and t to DOY and TOD
%   Detailed explanation goes here

yr0 = ST0(1);
mon0=ST0(2);
day0=ST0(3);
hr0=ST0(4);
min0=ST0(5);
sec0 = ST0(6);

yrpassed = floor(t/31536000);
yr = yr0 + yrpassed;
t = t-floor(t/31536000);

daypassed = floor(t/86400);
t = t-floor(t/86400);


febflag = abs(sign(mod(yr,4))-1);

switch mon0
    case 1
        doy0 = day0;
    case 2
        doy0 = day0 + 31;
    case 3
        doy0 = day0 + febflag + 59;
    case 4
        doy0 = day0 + febflag + 90;
    case 5
        doy0 = day0 + febflag + 120;
    case 6
        doy0 = day0 + febflag + 151;
    case 7
        doy0 = day0 + febflag + 181;
    case 8
        doy0 = day0 + febflag + 212;
    case 9
        doy0 = day0 + febflag + 243;
    case 10
        doy0 = day0 + febflag + 273;
    case 11
        doy0 = day0 + febflag + 304;
    case 12
        doy0 = day0 + febflag + 334;
    otherwise
        disp('trololo')
end

doy = doy0 + daypassed;
doy = mod(doy,365);


tod0 = hr0*3600 + min0*60 + sec0;
tod = tod0 + t;
tod = mod(tod,86400);

end

