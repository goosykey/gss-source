function [ F81,F107, Kp ] = getsolar( jd, sigmalevel )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[m,n]=size(jd);

Fb = 125; %mean mean mean F81
Fm = 40; %F81 variance amplitude
dS = 2451545 + 13.5*365; %shift to reference day from J2000.0
pS = pi/2; %sin phase shift (+pi/2 if max @ref.day, -pi/2 if min, 0 if grows)
T = 11.1; %years cycle
Kpmean = 3;

F81  = Fb + Fm*sin(2*pi*(jd-dS)/T/365 + pS);
F107 = F81 + 4.*sigmalevel.*randn(m,n);
Kp = round(Kpmean+randn(m,n));
if Kp<0
    Kp=0;
end
    
% 
% F81 = 125;
% F107 = 125;
% Kp = 3;


end

