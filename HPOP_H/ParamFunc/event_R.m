function [value, isterminal, direction] = event_R(T, Y)
%EVENT_SMA Summary of this function goes here
%   Detailed explanation goes here


R0 = 6378.14;

r = norm(Y(1:3));

value = (r < (R0+150));

isterminal = 1;

direction = 0;

end