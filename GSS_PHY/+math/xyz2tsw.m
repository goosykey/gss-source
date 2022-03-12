function [ Rtsw ] = xyz2tsw( Rxyz, O, i, u )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


M3 = [cos(pi/2-O) -sin(pi/2-O) 0
      sin(pi/2-O)  cos(pi/2-O) 0
        0               0      1];
    
M2 = [cos(pi-i) 0 sin(pi-i)
        0       1     0
     -sin(pi-i) 0 cos(pi-i)];
 
M1 = [cos(u) -sin(u) 0
      sin(u)  cos(u) 0
        0       0    1];
    
M = M1*M2*M3;

Rtsw = M*Rxyz;




end

