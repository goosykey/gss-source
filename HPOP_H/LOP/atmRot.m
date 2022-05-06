function [wa] = atmRot(H)
%ATMROT Summary of this function goes here
%   Detailed explanation goes here

w0 = 7.2921159e-05; %rad/s

h = [200 240 280 320 360 400 500];
L = [1.1 1.2 1.3 1.38 1.42 1.45 1.52];

hq = 150:800;

LL = interp1(h,L,hq, 'spline', 'extrap');

if H > 150
    wa = LL(round(H)==hq) * w0;
else
    wa = w0;
end

end