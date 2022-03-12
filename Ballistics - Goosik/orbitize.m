function [ Y ] = orbitize( y )
%CONVERTER a -> p
%   Converts a - e - om - Om - i - u to p - l1 - l2 - Om - i - u

Y = y;
Y(1) = y(1) * (1-y(2)^2);
Y(2) = y(2) * cos (y(3));
Y(3) = y(2) * sin (y(3));

end

