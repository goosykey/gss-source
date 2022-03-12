function [ a,e,om,Om,in,u ] = xyz2kepler1( x,y,z,vx,vy,vz )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

mu = 398600.4418; % Earth grav.param.

rxyz = [x,y,z]';
r = norm(rxyz);
vxyz = [vx,vy,vz]';

hxyz = cross(rxyz,vxyz); 
h = norm(hxyz);

exyz = 1/mu * (cross(vxyz,hxyz))-rxyz/r;
e = norm(exyz);

a = h^2/mu/(1-e^2);

kxyz = [0;0;1];
in = acos(dot(kxyz,hxyz)/h);

nxyz = cross(kxyz,hxyz);

Om = acos(dot([1;0;0],nxyz)/norm(nxyz));
if dot (nxyz,[0;1;0]) < 0
    Om = 2*pi-Om;
end

om = acos(dot(nxyz,exyz)/norm(nxyz)/e);
if dot (exyz,[0;0;1]) < 0
    om = 2*pi-om;
end

u = atan2(z/sin(in), x*cos(Om)+y*sin(Om));



end

