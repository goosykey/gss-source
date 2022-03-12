function [ inc ] = sunsyn( alt )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

mu    = 398600.440;      % Earth’s gravitational parameter [km^3/s^2]
Re = 6378;               % Earth radius [km]
J2  = 0.0010826269;      % Second zonal gravity harmonic of the Earth
we = 1.99106e-7;    % Mean motion of the Earth in its orbit around the Sun [rad/s]

a = alt + Re;
e=0;
h = a;
n = (mu./a.^3).^0.5;

tol = 1e-10;         % Error tolerance
% Initial guess for the orbital inclination
i0 = 180/pi*acos(-2/3*(h/Re).^2*we./(n*J2));
err = 1e1;
while(err >= tol )
    % J2 perturbed mean motion
    np  = n.*(1 + 1.5*J2*(Re./h).^2.*(1 - e^2)^0.5.*(1 - 3/2*sind(i0).^2));
    i = 180/pi*acos(-2/3*(h/Re).^2*we./(np*J2));
    err = abs(i - i0);
    i0 = i;
end

inc = i;


end

