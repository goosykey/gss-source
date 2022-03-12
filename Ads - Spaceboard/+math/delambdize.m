function [ ti,ephi ] = delambdize( L5i,yi )
%DELAMBDIZE Summary of this function goes here
%   Detailed explanation goes here

ti = yi(:,6);
ephi = yi;

ephi(:,1) = yi(:,1);
ephi(:,2) = sqrt(yi(:,2).^2 + yi(:,3).^2);
ephi(:,3) = atan2(yi(:,3),yi(:,2)) - atan2(yi(:,5),yi(:,4));
ephi(:,4) = atan2(yi(:,5),yi(:,4));
ephi(:,5) = 2*asin( sqrt(yi(:,4).^2 + yi(:,5).^2) );
ephi(:,6) = L5i - ephi(:,4);

end

