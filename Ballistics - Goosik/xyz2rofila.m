function [ Rrofila ] = xyz2rofila( Rxyz )
%ECEF -> rho, lat, lon 
%   Converts ECEF rotating (!) XYZ to ro, fi, la (spherical)

x = Rxyz(1);
y = Rxyz(2);
z = Rxyz(3);
ro = sqrt(x^2 + y^2 + z^2);
fi = asin(z/sqrt(x^2 + y^2 + z^2));

la = atan2(y,x);

Rrofila = [ro;fi;la];


end

