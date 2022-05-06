function Dy = earth353(t, y, varargin)
%>>>MAIN MOTION MODEL <<<
%   defines perturbations and motion equations
%
%
%     INPUTS:
%     t : nominal orbit parameters
%     y : single ephemeris point
%       
% 
%     OUTPUTS:
%     Dy : ephemeris DERIVATIVES
% 


%% VARARGIN DEFAULTS

jd0 = astro.jday(2020, 3, 20, 12, 0, 0);
f_Sm = @(jd, y)1; % default value of 1 m^2
f_Se = @(jd, y)1; % default value of 1 m^2
thruster = nan; % by default, thruster does not exist to save time creating objects
f_thrustflag = @(jd, y)0;
f_thrustVector = @(jd, y)[1,0,0];
f_Cxa = @(jd, y)2.2;


%% PARSE VARARGIN

while ~isempty(varargin)
    switch lower(varargin{1})
        case 'jd0'
            jd0 = varargin{2};
        case 'f_sm'
            f_Sm = varargin{2};
        case 'f_se'
            f_Se = varargin{2};
        case 'f_thrustflag'
            f_thrustflag = varargin{2};
        case 'f_thrustvector'
            f_thrustVector = varargin{2};
        case 'f_cxa'
            f_Cxa = varargin{2};
        otherwise
            %varargin{1}
            error(['Unexpected option: ' varargin{1}])
    end
    varargin(1:2) = [];
end

%% INITIALIZE

Dy = y(:);

p = y(1);
l1 = y(2);
l2 = y(3);
Om = y(4);
in = y(5);
u = y(6);
m = y(7);

jd = jd0 + t/86400;
[month, day, ~] = astro.gdate (jd);
DOY = 30*(month-1) + round(day);


%% PERTURBATIONS ON/OFF - PRECISION OPTIONS
       
PREC = 2;
switch PREC
    case 0
        Zg = 0;         % Advanced gravity nonspherics
        Zj2 = 1;        % Simplified (J2 only) nonspherics
        Zatm = 0;       % Atmospheric drag
        Zsunmoon = 0;   % Sun & Moon gravity perturbation
        Zemp = 0;       % Electromagnetic pressure (Sun)
        Zprop = 1;      % Propulsion
        maxdeg = 0;     % Max. degree for nonspherics
        maxord = 0;     % Max. order for nonspherics
        precnut = 0;    % Earth axis precession & nutation for nonspherics
    case 1
        Zg = 1;
        Zj2 = 0;
        Zatm = 1;
        Zsunmoon = 1;
        Zemp = 0;
        Zprop = 1;        
        maxdeg = 4;
        maxord = 0;
        precnut = 0;
    case 2
        Zg = 1;
        Zj2 = 0;
        Zatm = 0;
        Zsunmoon = 0;
        Zemp = 0;
        Zprop = 0;        
        maxdeg = 2;
        maxord = 2;
        precnut = 1;       
end

%% CONSTANTS
del = 66.07354399E03; % includes J2
mu = 398600.4415; % Earth grav.param. km3/s2
muM = 7.36E+22 * 6.67384E-20; % Moon grav.param. km3/s2
muS = 1.98892E+30 * 6.67384E-20; % Sun grav.param. km3/s2
au1 = 149597871; % 1 Astronomic Unit in kms
R0 = 6378.14; % mean Earth radius
g0 = 9.80665 * 10^(-3); % g0 km/s2
I0 = 1366; %W/m^2
c = 3E+08; %m/s

%% BASIC CALCULATIONS

%RADIUS
e = sqrt(l1.^2 + l2.^2);
om = atan2(l2,l1);
nu = u-om;
r = p / (1 + e*cos(nu));
a = p/(1-e^2);

%Orbital rvector thus is [0, r, 0]'
%rxyz = math.tsw2xyz([0, r, 0]',Om,in,u); %rxyz J2000


%% R and V

[rxyz, vxyz] = math.randv(a,e,in,Om,om,nu);
v = sqrt(sum(vxyz));


%% SUNMOON

[rsun,rtascsun,declsun] = astro.sun ( jd ); %Vallado

if Zsunmoon ~= 0
        
    r12 = astro.moon(jd)'*R0; r12abs = sqrt(sum(r12 .* r12));
    r13 = rsun*au1; r13abs = sqrt(sum(r13 .* r13));
    r2 = rxyz - r12;    r2abs = sqrt(sum(r2 .* r2));
    r3 = rxyz - r13;    r3abs = sqrt(sum(r3 .* r3));
    
    %XYZ cartesian 3body perturbations
    fsunxyz =  -muS/r3abs^3*r3 - muS/r13abs^3*r13; fsuntsw  = math.xyz2tsw(fsunxyz, Om,in,u);
    fmoonxyz = -muM/r2abs^3*r2 - muM/r12abs^3*r12; fmoontsw = math.xyz2tsw(fmoonxyz,Om,in,u);
    Tsunmoon = fsuntsw(1) + fmoontsw(1);
    Ssunmoon = fsuntsw(2) + fmoontsw(2);
    Wsunmoon = fsuntsw(3) + fmoontsw(3);
else
    Tsunmoon = 0;
    Ssunmoon = 0;
    Wsunmoon = 0;    
end

%% DRAG & EMP PARAMETERS

Sm = f_Sm(jd, y);
Semp = f_Se(jd, y);
Cxa = f_Cxa(jd, y);

%% EM PRESSURE

los = astro.sight(rxyz,rsun*au1,'e');
if strcmp(los,'no ')
    Zemp = 0;
end

if Zemp ~=0
    pem = I0/c*2;
    fempxyz = -pem*Semp/m*rsun/norm(rsun')/10^3; %<<km!
    femptsw = math.xyz2tsw(fempxyz,Om,in,u);
    Temp = femptsw(1);
    Semp = femptsw(2);
    Wemp = femptsw(3);
else 
    Temp = 0;
    Semp = 0;
    Wemp = 0;
end




%% nonspherics full

if precnut ~= 0
    ttt= (jd - 2451545.0  )/ 36525.0; %Julian Cent
    rteme = astro.eci2teme(rxyz,[0;0;0],[0;0;0],ttt,106,0,'a'); %rvec J2000 to TEME
else
    rteme = rxyz; %no prec nut
end

recef = math.R3(astro.gstime(jd))*rteme;

if Zg ~= 0
    fgxyz = grav (recef, jd, maxdeg, maxord, precnut, false);
    fgtsw = math.xyz2tsw(fgxyz,Om,in,u);
    Tg = fgtsw(1);
    Sg = fgtsw(2);
    Wg = fgtsw(3);
else
    Tg = 0;
    Sg = 0;
    Wg = 0;
end



%% nonspherics simplified (J2 Perturbation)
if Zj2 ~= 0 && Zg ==0
    Tj2 = -mu*del/r^4 * sin(2*u) * (sin(in))^2;
    Sj2 = mu*del/r^4 * (3*(sin(u))^2*(sin(in))^2-1);
    Wj2 = mu*del/r^4 * sin(u) * sin(2*in);
else
    Tj2 = 0;
    Sj2 = 0;
    Wj2 = 0;
end

%% ATMOSPHERIC DRAG
if Zatm ~=0
    [ F81,F107, Kp ] = astro.getsolar( jd, 3 );
    rho = astro.atmgost2(recef, r-6378.14,DOY,astro.gstime(jd),rtascsun,declsun,F107,F81,Kp); % GOST
    %rho = astro.atmgost2(recef, r-6378.14,DOY,astro.gstime(jd),rtascsun,declsun,125,125,2); % GOST
    Tatm = -Cxa*rho*(v*1000)^2*Sm/(2*m) * 10^(-3);
else
    Tatm = 0;
end

%% PROPULSION

if ~isnan(thruster)
    thrust = thruster.thrust;
    Isp = thruster.Isp;
else
    thrust = 0;
    Isp = 1550;
end

thrustVector = f_thrustVector(jd, y);
thrustflag = f_thrustflag(jd, y);

thrust = thrust * 10^(-3); % m to km
fprop = thrust/m;
mdot = thrust/(Isp*g0);

fprop_xyz = fprop*thrustVector;

Tprop = fprop_xyz(1);
Sprop = fprop_xyz(2);
Wprop = fprop_xyz(3);

%% TOTAL PERTURBATIONS
T = Zj2*Tj2 + Zg*Tg + Zsunmoon*Tsunmoon + Zemp*Temp + Zprop*thrustflag*Tprop + Zatm*Tatm;
S = Zj2*Sj2 + Zg*Sg + Zsunmoon*Ssunmoon + Zemp*Semp + Zprop*thrustflag*Sprop;
W = Zj2*Wj2 + Zg*Wg + Zsunmoon*Wsunmoon + Zemp*Wemp + Zprop*thrustflag*Wprop;


%% FINALLY -- EQUATIONS
Dy(1) = 2*r*(p/mu)^(1/2)*T;
Dy(2) = sqrt(p/mu)*(T*(1+r/p)*cos(u) + S*sin(u) + r/p*(T*l2-W*l1/tan(in)*sin(u)));
Dy(3) = sqrt(p/mu)*(T*(1+r/p)*sin(u) - S*cos(u) + r/p*(T*l1+W*l2/tan(in)*sin(u)));
Dy(4) = -r/(mu*p)^(1/2)*sin(u)/sin(in)*W;
Dy(5) = -r*cos(u)/(mu*p)^(1/2)*W;
Dy(6) = (mu*p)^(1/2)/r^2 + r/(mu*p)^(1/2)*sin(u)/tan(in)*W;
Dy(7) = -Zprop*thrustflag*mdot;   


end


function fgxyz = grav (recef, jd, maxdeg, maxord, precnut, useAerospace)

ttt= (jd - 2451545.0  )/ 36525.0; %Julian Cent
R0 = 6378.137;
mu = 398600.4415;

if useAerospace
    norm(recef'*1e3)
    [gx, gy, gz] = gravitysphericalharmonic(recef'*1e3,'EGM2008',2);
    fgecef = [gx, gy, gz]'/1e3;
else
    rrofila = math.xyz2rofila(recef); %ECEF to Rho, lat, lon
    RO = rrofila(1); dUro = 0;
    FI = rrofila(2); dUfi = 0;
    LA = rrofila(3); dUla = 0;
    
    LEG = zeros(maxdeg-1,maxdeg+2);
    for q = 2:maxdeg
        LEG(q-1,1:q+1)=(legendre(q,sin(FI)))';
    end
    
    mtP  = math.mtanleg( maxdeg,maxdeg,FI );
    
    for q = 2:maxdeg
        for w = 0:min([q,maxord])
            fuck = (-1)^w; %KONDON FACTOR QUANTUM CRAP EXCLUSION
            uro = (R0/RO)^q*(q+1)*fuck*LEG(q-1,w+1)*(astro.gcons('C',q,w)*cos(w*LA)+astro.gcons('S',q,w)*sin(w*LA));
            ufi = (R0/RO)^q*(fuck*(-1)*LEG(q-1,w+2)-mtP(q+1,w+1))*(astro.gcons('C',q,w)*cos(w*LA)+astro.gcons('S',q,w)*sin(w*LA));
            ula = (R0/RO)^q*w*fuck*LEG(q-1,w+1)*(astro.gcons('S',q,w)*cos(w*LA)-astro.gcons('C',q,w)*sin(w*LA));
            
            dUro = dUro - mu/RO^2*uro;
            dUfi = dUfi + mu/RO*ufi;
            dUla = dUla + mu/RO*ula;
        end
    end
    
    fgecef = [(dUro/RO - 1/sqrt(recef(1)^2+recef(2)^2)*recef(3)/RO^2*dUfi)*recef(1) - dUla * recef(2)/(recef(1)^2+recef(2)^2);
        (dUro/RO - 1/sqrt(recef(1)^2+recef(2)^2)*recef(3)/RO^2*dUfi)*recef(2) + dUla * recef(1)/(recef(1)^2+recef(2)^2);
        1/RO*dUro*recef(3) + sqrt(recef(1)^2+recef(2)^2)/RO^2*dUfi];
    
end


fgteme = math.R3(-astro.gstime(jd))*fgecef;

if precnut ~= 0
    fgxyz = astro.teme2eci(fgteme,[0;0;0],[0;0;0],ttt,106,0,'a');
else
    fgxyz = fgteme; %no prec nut
end



end