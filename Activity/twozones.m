function [ Rurb, Rrur ] = twozones( pop )
%TWOZONES Summary of this function goes here
%   Detailed explanation goes here

p1 = [-0.6667   12.8333  -87.3333  253.1667 -266.0000];
p2 = [-0.5625   12.3750  -93.1875  296.3750 -337.5000];

a = polyval(p1, log10(pop));
b = polyval(p2, log10(pop));

Rurb = (abs(a)+a)/2 * 1.2 + 1;
Rrur = (abs(b)+b)/2 * 1.3 + 1;



end

