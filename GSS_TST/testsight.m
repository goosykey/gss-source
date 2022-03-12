function [d] = testsight(r, rs)
%TESTSIGHT Summary of this function goes here
%   Detailed explanation goes here

R0 = 6378.14;

D = dot(r, rs, 2).^2 - (dot(r,r,2) - R0.^2);

d1 = - dot(r, rs, 2) + sqrt(D);

d = D < 0 | d1 < 0;

end

