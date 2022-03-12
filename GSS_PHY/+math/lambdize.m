function [ Y ] = lambdize( y )
%CONVERTER a -> p
%   Converts a - e - om - Om - i - u to L0 L1 L2 L3 L4 L5

Y = y;
Y(1) = y(1);
Y(2) = y(2) * cos (y(3)+y(4));
Y(3) = y(2) * sin (y(3)+y(4));
Y(4) = sin(y(5)/2) * cos(y(4));
Y(5) = sin(y(5)/2) * sin(y(4));
Y(6) = y(6) + y(4);

Y(7) = y(7);

end

