function [ rtopo ] = xyz2topo( rxyz, fignd, lagnd, jd )
%J2000 XYZ to GND STATION XYZ
%   Detailed explanation goes here

st = gstime(jd) + lagnd;
rho = 6378.14;


rgndxyz = rho*[cos(fignd)*cos(st);cos(fignd)*sin(st);sin(fignd)];

e1 = rgndxyz/norm(rgndxyz);
e3 = [-sin(fignd)*cos(st);-sin(fignd)*sin(st);cos(fignd)];
e2 = cross(e3,e1);

B = [e1';e2';e3'];
C = B';

%rtopo = C\(rxyz)-rgndxyz;
rtopo = C\(rxyz-rgndxyz);


end

