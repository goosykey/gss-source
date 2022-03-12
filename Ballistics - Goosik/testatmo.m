function [ MAM ] = testatmo( N )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


R0 = 6378.14; % mean Earth radius
mu = 398600.44; % Earth grav.param. km3/s2

recef = zeros(N,3);
jd = zeros(N,1);
gst = zeros(N,1);
DOY = zeros(N,1);
h = zeros(N,1);
alph = zeros(N,1);
delt = zeros(N,1);
solar = zeros(N,3);
KK = zeros(N,5);
RhoN = zeros(N,1);
Rho = zeros(N,1);
Rhofu = zeros(N,1);

a = 3000;
b = 6000;



for i = 1:N
    clc
    fprintf('%3.2f%% COMPLETE;\n', (i-1)/N*100);
    sig = sign(rand-0.5);
    while (norm(recef(i,:)) - 6378.14) < 120 || (norm(recef(i,:)) - 6378.14) >1500
        recef(i,:) = (a + (b-a)*rand(1,3))*sig;        
    end
    jd(i) = 2457389 + 365 * rand;
    if i == 1
        recef(i,:) = [-5586.80663515860 -3038.28702598473 -4196.54126151769];
        jd(i) = [2457631.44646356];        
    end
    gst(i) = gstime(jd(i));
    DOY(i) = rem(jd(i)-2457388.5,365);
    h(i)=norm(recef(i,:))-6378.137;
    [rsun,alph(i),delt(i)]=sun(jd(i));
    [solar(i,1), solar(i,2), solar(i,3)] = getsolar(jd(i),3);
    if i == 1
        solar(i,1) = 98.1882597744764;
        solar(i,2) = 111.016023748347;
        solar(i,3) = 5;
    end
    [Rho(i),RhoN(i),KK(i,:)]=atmgost2(recef(i,:),h(i),DOY(i),gst(i),alph(i),delt(i),solar(i,2), solar(i,1), solar(i,3));
    Rhofu(i) = atmgost(recef(i,:),h(i),DOY(i),0,gst(i),alph(i),delt(i),solar(i,2), solar(i,1), solar(i,3));
    
    
end

clc;
fprintf('TRANSFORMATION COMPLETE;\n');

MAM = [jd,recef,h,DOY,gst,alph,delt,solar,RhoN,KK,Rho,Rhofu];

end

