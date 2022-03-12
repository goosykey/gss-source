function [ dUro, dUfi, dUla ] = rofila1test( RO, FI, LA )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

R0 = 6378.137; % mean Earth radius
mu = 398600.4418;

dUro = 0;
dUfi = 0;
dUla = 0;

LEG = zeros(7,10);
for q = 2:8
    LEG(q-1,1:q+1)=(legendre(q,sin(FI)))';
end

mtP  = mtanleg( 8,8,FI );

for q = 2:8
    for w = 0:q
        fuck = (-1)^w; %KONDON FACTOR QUANTUM CRAP EXCLUSION
        uro = (R0/RO)^q*(q+1)*fuck*LEG(q-1,w+1)*(gcons('C',q,w)*cos(w*LA)+gcons('S',q,w)*sin(w*LA));
        ufi = (R0/RO)^q*(fuck*(-1)*LEG(q-1,w+2)-mtP(q+1,w+1))*(gcons('C',q,w)*cos(w*LA)+gcons('S',q,w)*sin(w*LA));
        ula = (R0/RO)^q*w*fuck*LEG(q-1,w+1)*(gcons('S',q,w)*cos(w*LA)-gcons('C',q,w)*sin(w*LA));
        
        dUro = dUro - mu/RO^2*uro;
        dUfi = dUfi + mu/RO*ufi;
        dUla = dUla + mu/RO*ula;
    end
end

end

