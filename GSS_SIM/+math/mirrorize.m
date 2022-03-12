function [ CONNM ] = mirrorize( CONN )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

CONNM = CONN + CONN';

for i = 1 : length(CONNM(1,:))
    CONNM(i,i) = CONNM(i,i)/2;
end


end

