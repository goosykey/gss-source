function [value, isterminal, direction] = event_SMA(T, Y)
%EVENT_SMA Summary of this function goes here
%   Detailed explanation goes here

@(a)T;

R0 = 6378.14;

p = Y(1);
l1 = Y(2);
l2 = Y(3);
u = Y(6);

e = sqrt(l1.^2 + l2.^2);
om = atan2(l2,l1);
nu = u-om;
r = p / (1 + e*cos(nu));

value = (r < (R0+150));

isterminal = 1;

direction = 0;

end