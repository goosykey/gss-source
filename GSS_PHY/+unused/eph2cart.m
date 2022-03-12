function [ EPHXYZ ] = eph2cart( EPH )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

yc = unorbitize(EPH);

a = yc (:,1);
e = yc (:,2);
om = yc(:,3);
Om = yc(:,4);
in = yc(:,5);
u = yc (:,6);
m = yc (:,7);

nu = u-om;

[R,V] = randv(a,e,in,Om,om,nu);

EPHXYZ = [R',V',m];



end