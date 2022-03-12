function [ dUxyz ] = rofila2xyzJAC( dUrofila, ro, fi, la, jd )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

thetag = gstime(jd);

JAC = [cos(fi)*cos(la+thetag) -ro*sin(fi)*cos(la+thetag) -ro*cos(fi)*sin(la+thetag)
    cos(fi)*sin(la+thetag) -ro*sin(fi)*sin(la+thetag) ro*cos(fi)*cos(la+thetag)
    sin(fi) cos(fi) 0];

dUxyz = JAC'\dUrofila;    


end

