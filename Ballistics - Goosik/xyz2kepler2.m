function [ a,e,om,Om,in,u ] = xyz2kepler2( x,y,z,vx,vy,vz )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

mu = 398600.4418; % Earth grav.param.

R = [x,y,z]; r = norm(R);
V = [vx,vy,vz]; v = norm(V);

H = cross(R,V); h = norm(H);

Om = atan2(H(1),-H(2));
in = atan2(sqrt(H(1)^2+H(2)^2),H(3));

Om = mod(Om,2*pi);
in = mod(in,2*pi);

P = R1(in)*R3(Om)*R';
u = atan2(P(2),P(1));
u = mod(u,2*pi);

a = mu*r/(2*mu-r*v^2);
e = sqrt(1-h^2/mu/a);

cosE = (a-r)/a/e;
sinE = dot(R,V)/e/sqrt(mu*a);

nu = atan2(sqrt(1-e^2)*sinE,cosE-e);

om = u-nu;
om = mod(om,2*pi);



end

