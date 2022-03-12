function [ V1 ] = orderize( V )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

V1=intersect(V,V);
V1 = V1(V1~=0);



end

