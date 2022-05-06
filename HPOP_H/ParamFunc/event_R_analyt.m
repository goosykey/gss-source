function [value, isterminal, direction] = event_R_analyt(T, y)
%EVENT_SMA Summary of this function goes here
%   Detailed explanation goes here

a = y(1);
e = y(2);
M = y(6);

E = M;
delta = inf;

while delta > 1e-06
    E1 = M + e.*sin(E);
    delta = abs(E-E1);
    E = E1;
end

nu = atan2( sqrt(1-e.^2).*sin(E),(cos(E)-e) );
p = a*(1-e*e);
r = p / (1 + e*cos(nu));

R0 = 6378.14;

value = (r < (R0+150));

isterminal = 1;

direction = 0;

end
