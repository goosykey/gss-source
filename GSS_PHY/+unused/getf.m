function [ fT, fW ] = getf( XSI, u, m )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

ulim = 55 * pi/180;
P = 0.179 * 10^(-3);


if abs(sin(u)) < abs(sin(ulim))
    fT = cos(XSI)*P/m;
    fW = -sign(cos(u))*sin(XSI)*P/m;
else
    fT = 0;
    fW = 0;
end

end

