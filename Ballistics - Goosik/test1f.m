function [ fgtsw, fatmtsw, fsunmoontsw, femptsw, ftsw ] = test1f( p, l1, l2, Om, in, u, m, jd, F107, F81, Kp )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

ttt= (jd - 2451545.0  )/ 36525.0; %Julian Cent
DOY = rem(jd-2457388.5,365);

maxdeg = 8;
maxord = 8;

%CONSTANTS
mu = 398600.44; % Earth grav.param. km3/s2
muM = 7.36E+22 * 6.67384E-20; % Moon grav.param. km3/s2
muS = 1.98892E+30 * 6.67384E-20; % Sun grav.param. km3/s2
au1 = 149597871; % 1 Astronomic Unit in kms
R0 = 6378.14; % mean Earth radius
I0 = 1366; %W/m^2
c = 3E+08; %m/s

e = sqrt(l1.^2 + l2.^2);
om = perigee(l1,l2);
nu = u-om;
r = p / (1 + e*cos(nu));

Cxa = 2.2;
Sm = 20;
Semp = 20;

%Orbital rvector thus is [0, r, 0]'
rxyz = tsw2xyz([0, r, 0]',Om,in,u); %rxyz J2000
[rsun,rtascsun,declsun] = sun ( jd ); %Vallado


%SUNMOON
r12 = moon(jd)'*R0; r12abs = sqrt(sum(r12 .* r12));
r13 = rsun'*au1; r13abs = sqrt(sum(r13 .* r13));
r2 = rxyz - r12;    r2abs = sqrt(sum(r2 .* r2));
r3 = rxyz - r13;    r3abs = sqrt(sum(r3 .* r3));
%XYZ cartesian 3body perturbations
fsunmoonxyz = -muM/r2abs^3*r2 - muM/r12abs^3*r12 ...
    - muS/r3abs^3*r3 - muS/r13abs^3*r13;
fsunmoontsw = xyz2tsw(fsunmoonxyz,Om,in,u);

%EMP
pem = I0/c*2;
fempxyz = -pem*Semp/m*rsun'/norm(rsun)/10^3; %<<km!
femptsw = xyz2tsw(fempxyz,Om,in,u);

%velocity
v = sqrt(mu*(2/r - (1-e^2)/p));

%g
rteme = eci2teme(rxyz,[0;0;0],[0;0;0],ttt,106,0,'a'); %rvec J2000 to TEME
recef = R3(gstime(jd))*rteme;

rrofila = xyz2rofila(recef); %ECEF to Rho, lat, lon
RO = rrofila(1); dUro = 0;
FI = rrofila(2); dUfi = 0;
LA = rrofila(3); dUla = 0;

LEG = zeros(maxdeg-1,maxdeg+2);
for q = 2:maxdeg
    LEG(q-1,1:q+1)=(legendre(q,sin(FI)))';
end

mtP  = mtanleg( maxdeg,maxdeg,FI );

for q = 2:maxdeg
    for w = 0:min([q,maxord])
        fuck = (-1)^w; %KONDON FACTOR QUANTUM CRAP EXCLUSION
        uro = (R0/RO)^q*(q+1)*fuck*LEG(q-1,w+1)*(gcons('C',q,w)*cos(w*LA)+gcons('S',q,w)*sin(w*LA));
        ufi = (R0/RO)^q*(fuck*(-1)*LEG(q-1,w+2)-mtP(q+1,w+1))*(gcons('C',q,w)*cos(w*LA)+gcons('S',q,w)*sin(w*LA));
        ula = (R0/RO)^q*w*fuck*LEG(q-1,w+1)*(gcons('S',q,w)*cos(w*LA)-gcons('C',q,w)*sin(w*LA));
        
        dUro = dUro - mu/RO^2*uro;
        dUfi = dUfi + mu/RO*ufi;
        dUla = dUla + mu/RO*ula;
    end
end

fgecef = [(dUro/RO - 1/sqrt(recef(1)^2+recef(2)^2)*recef(3)/RO^2*dUfi)*recef(1) - dUla * recef(2)/(recef(1)^2+recef(2)^2);
    (dUro/RO - 1/sqrt(recef(1)^2+recef(2)^2)*recef(3)/RO^2*dUfi)*recef(2) + dUla * recef(1)/(recef(1)^2+recef(2)^2);
    1/RO*dUro*recef(3) + sqrt(recef(1)^2+recef(2)^2)/RO^2*dUfi];

fgteme = R3(-gstime(jd))*fgecef;
fgxyz = teme2eci(fgteme,[0;0;0],[0;0;0],ttt,106,0,'a');
fgtsw = xyz2tsw(fgxyz,Om,in,u);

%ATM
rho = atmgost2(recef, r-6378.14,DOY,gstime(jd),rtascsun,declsun,F107,F81,Kp); % GOST
Tatm = -Cxa*rho*(v*1000)^2*Sm/(2*m) * 10^(-3);
fatmtsw = [Tatm;0;0];

ftsw = fgtsw + fatmtsw + femptsw + fsunmoontsw;



end

