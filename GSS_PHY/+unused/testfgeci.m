function [MMM] = testfgeci( N )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

R0 = 6378.14; % mean Earth radius
mu = 398600.44; % Earth grav.param. km3/s2

reci = zeros(N,3);
rteme = zeros(N,3);
recef = zeros(N,3);
rrofila = zeros(N,3);
jd = zeros(N,1);
gst = zeros(N,1);
fgeci = zeros(N,3);
fgteme = zeros(N,3);
fgecef = zeros(N,3);
dUs = zeros(N,3);

a = 3000;
b = 6000;



for i = 1:N
    clc
    fprintf('%3.2f%% COMPLETE;\n', (i-1)/N*100);
    sig = sign(rand-0.5);
    reci(i,:) = (a + (b-a)*rand(1,3))*sig;
    jd(i) = 2451545.0 + 12000 * rand;
    if i == 1
        reci(i,:) = [5780.57872015653, 3024.88697217792, 5473.88502799918];
        jd(i) = 2460753.03041633;        
    end
    gst(i) = gstime(jd(i));
    ttt = (jd(i) - 2451545.0  )/ 36525.0; %Julian Cent
    rteme(i,:) = eci2teme(reci(i,:)',[0;0;0],[0;0;0],ttt,106,0,'a');
    
    recef(i,:) = R3(gstime(jd(i)))*rteme(i,:)';
    
    rrofila(i,:) = xyz2rofila(recef(i,:)); %ECEF to Rho, lat, lon
    RO = rrofila(i,1); dUro = 0;
    FI = rrofila(i,2); dUfi = 0;
    LA = rrofila(i,3); dUla = 0;
    
    
     LEG = zeros(7,10);
    for q = 2:8
        LEG(q-1,1:q+1)=(legendre(q,sin(FI)))';
    end
    
    mtP  = mtanleg( 8,8,FI );
    
    for q = 2:8
        for w = 0:8
            fuck = (-1)^w; %KONDON FACTOR QUANTUM CRAP EXCLUSION
            uro = (R0/RO)^q*(q+1)*fuck*LEG(q-1,w+1)*(gcons('C',q,w)*cos(w*LA)+gcons('S',q,w)*sin(w*LA));
            ufi = (R0/RO)^q*(fuck*(-1)*LEG(q-1,w+2)-mtP(q+1,w+1))*(gcons('C',q,w)*cos(w*LA)+gcons('S',q,w)*sin(w*LA));
            ula = (R0/RO)^q*w*fuck*LEG(q-1,w+1)*(gcons('S',q,w)*cos(w*LA)-gcons('C',q,w)*sin(w*LA));
            
            dUro = dUro - mu/RO^2*uro;
            dUfi = dUfi + mu/RO*ufi;
            dUla = dUla + mu/RO*ula;
            dUs (i,:) = [dUro,dUfi,dUla];
        end
    end
    
    fgecef(i,:) = [(dUro/RO - 1/sqrt(recef(i,1)^2+recef(i,2)^2)*recef(i,3)/RO^2*dUfi)*recef(i,1) - dUla * recef(i,2)/(recef(i,1)^2+recef(i,2)^2);
        (dUro/RO - 1/sqrt(recef(i,1)^2+recef(i,2)^2)*recef(i,3)/RO^2*dUfi)*recef(i,2) + dUla * recef(i,1)/(recef(i,1)^2+recef(i,2)^2);
        1/RO*dUro*recef(i,3) + sqrt(recef(i,1)^2+recef(i,2)^2)/RO^2*dUfi];
    
    fgteme(i,:) = R3(-gstime(jd(i)))*fgecef(i,:)';
    
    fgeci(i,:) = teme2eci(fgteme(i,:)',[0;0;0],[0;0;0],ttt,106,0,'a');
end

clc;
fprintf('TRANSFORMATION COMPLETE;\n');

MMM = [reci, jd, gst, rteme, recef, rrofila, dUs,fgecef,fgteme,fgeci];




end

