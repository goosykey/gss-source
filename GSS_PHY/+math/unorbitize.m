function [ Y ] = unorbitize( y )
%CONVERTER p -> a
%   Converts p - l1 - l2 - Om - i - u to a - e - om - Om - i - u

Y = y;

Y(:,2) = sqrt(y(:,2).^2 + y(:,3).^2);
Y(:,3) = atan2 (y(:,2),y(:,3));

Y(:,1) = y(:,1) ./ (1-Y(:,2).^2);


end

