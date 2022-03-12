function [ Rxyz ] = tsw2xyz( Rtsw, O, i, u )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%{
su = sin(u);
cu = cos(u);
sO = sin(O);
cO = cos(O);
si = sin(i);
ci = cos(i);
%}

%Матрица 
M3 = [cos(pi/2-O) -sin(pi/2-O) 0
      sin(pi/2-O)  cos(pi/2-O) 0
         0            0        1];
     
M2 = [cos(pi-i) 0 sin(pi-i)
           0    1     0
     -sin(pi-i) 0 cos(pi-i)];
 
M1 = [cos(u) -sin(u) 0
      sin(u) cos(u)  0
         0      0    1];
     
M = M1*M2*M3;

Rxyz = M\Rtsw;

%Матрица по книге, 3.41:
%{ 
Mkn = [-su*cO-cu*sO*ci cu*cO-su*sO*ci -si*sO
    -su*sO+cu*cO*ci cu*sO+su*cO*ci si*cO
    cu*si su*si -ci];
Rxyzkn = Mkn*Rtsw
%}


end

