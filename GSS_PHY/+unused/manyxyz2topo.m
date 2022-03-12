function [ Rtopo ] = manyxyz2topo( T, Rxyz, fi, la )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

L = length(T);
global STARTTIME;

ST0 = STARTTIME;
jd0 = jday(ST0(1), ST0(2), ST0(3), ST0(4), ST0(5), ST0(6));

Rtopo = Rxyz-Rxyz;

for i = 1:L
    jd = jd0 + T(i)/86400;
    Rtopo(:,i) = xyz2topo(Rxyz(:,i),fi,la,jd);
end

end

