function Dy = earth353_FULL(t, y)
%>>>MAIN MOTION MODEL <<<
%   defines perturbations and motion equations
%
%   by A.Kharlan, © Yaliny, 2015
%
%     INPUTS:
%     t : nominal orbit parameters
%     y : single ephemeris point
% 
%     OUTPUTS:
%     Dy : ephemeris DERIVATIVES
% 
%     GLOBAL VARIABLES USED:
%     PROPXSI   : yaw XSI-angle for thrust calculation, defined elsewhere
%     PREC      : defined elsewhere, this determines the modelling precision
%     STARTTIME : defined elsewhere, this determines start date and UTC
%     BOUND     : defined elsewhere, this is sin(ulim) to limit active
%               stage during each orbit. Negative value means phasing,
%               positive value means stationkeeping;

Dy = y(:);

p = y(1);
l1 = y(2);
l2 = y(3);
Om = y(4);
in = y(5);
u = y(6);
m = y(7);


global PROPXSI BOUND PREC STARTTIME SATNAME; %TESTR;

if (BOUND >= 0 && abs(sin(u)) > BOUND) || (BOUND < 0 && abs(sin(u)) < -BOUND) % 
    thrustflag = 1;
else
    thrustflag = 0;
end

% if sin(u) >= 0
%     if cos(u) >= 0
%         PROPXSI = pi/2;
%         thrustflag = 1;
%     else
%         thrustflag = 0;
%         PROPXSI = -1;
%     end
% else
%     if cos(u) < 0
%         PROPXSI = pi/2;
%         thrustflag = 1;
%     else
%         thrustflag = 0;
%         PROPXSI = -1;
%     end
% end

%START TIME
ST0 = STARTTIME;
jd0 = jday(ST0(1), ST0(2), ST0(3), ST0(4), ST0(5), ST0(6));
jd = jd0 + t/86400;
[DOY,~] = doytod(ST0,t);
ttt= (jd - 2451545.0  )/ 36525.0; %Julian Cent

%PERTURBATIONS ON/OFF - PRECISION OPTIONS
       

switch PREC
    case 0
        Zg = 0;         % Advanced gravity nonspherics
        Zj2 = 1;        % Simplified (J2 only) nonspherics
        Zatm = 1;       % Atmospheric drag
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
        Zatm = 1;
        Zsunmoon = 1;
        Zemp = 1;
        Zprop = 1;        
        maxdeg = 8;
        maxord = 8;
        precnut = 1;       
end

if PROPXSI == -1
    Zprop = 0;
end


%CONSTANTS
del = 66.07E03; % includes J2
mu = 398600.44; % Earth grav.param. km3/s2
muM = 7.36E+22 * 6.67384E-20; % Moon grav.param. km3/s2
muS = 1.98892E+30 * 6.67384E-20; % Sun grav.param. km3/s2
au1 = 149597871; % 1 Astronomic Unit in kms
R0 = 6378.14; % mean Earth radius
g0 = 9.80665 * 10^(-3); % g0 km/s2
I0 = 1366; %W/m^2
c = 3E+08; %m/s

%RADIUS
e = sqrt(l1.^2 + l2.^2);
om = perigee(l1,l2);
nu = u-om;
r = p / (1 + e*cos(nu));

%Orbital rvector thus is [0, r, 0]'
rxyz = tsw2xyz([0, r, 0]',Om,in,u); %rxyz J2000
[rsun,rtascsun,declsun] = sun ( jd ); %Vallado

%SUNMOON
if Zsunmoon ~= 0
        
    r12 = moon(jd)'*R0; r12abs = sqrt(sum(r12 .* r12));
    r13 = rsun'*au1; r13abs = sqrt(sum(r13 .* r13));
    r2 = rxyz - r12;    r2abs = sqrt(sum(r2 .* r2));
    r3 = rxyz - r13;    r3abs = sqrt(sum(r3 .* r3));
    
    %XYZ cartesian 3body perturbations
    fsunmoonxyz = -muM/r2abs^3*r2 - muM/r12abs^3*r12 ...
        - muS/r3abs^3*r3 - muS/r13abs^3*r13;
    fsunmoontsw = xyz2tsw(fsunmoonxyz,Om,in,u);
    Tsunmoon = fsunmoontsw(1);
    Ssunmoon = fsunmoontsw(2);
    Wsunmoon = fsunmoontsw(3);
else
    Tsunmoon = 0;
    Ssunmoon = 0;
    Wsunmoon = 0;    
end

%DRAG & EMP PARAMETERS
if strcmp (SATNAME, 'beacon')
    Cxa = 1.1;%* 0.5*(sin(t/60)+1);
else
    Cxa = 2.2;
end

if strcmp (SATNAME, 'beacon')
    Sm = 3.9;
else
    nxyz = rsun'/norm(rsun');
    ntsw = xyz2tsw(nxyz,Om,in,u);
    cospsi = ntsw' * [1;0;0];
    Sm = 1.68 + 36.96 * abs(cospsi);
    Sm = 0.8;
end

if strcmp (SATNAME, 'beacon')
    Semp = 3.9;
else
    Semp = 36;
    Semp = 0.8;
end

%EM PRESSURE
if Zemp ~=0
    pem = I0/c*2;
    fempxyz = -pem*Semp/m*rsun'/norm(rsun)/10^3; %<<km!
    femptsw = xyz2tsw(fempxyz,Om,in,u);
    Temp = femptsw(1);
    Semp = femptsw(2);
    Wemp = femptsw(3);
else 
    Temp = 0;
    Semp = 0;
    Wemp = 0;
end



%velocity
v = sqrt(mu*(2/r - (1-e^2)/p));


%nonspherics full
if Zg ~= 0
    if precnut ~= 0
        rteme = eci2teme(rxyz,[0;0;0],[0;0;0],ttt,106,0,'a'); %rvec J2000 to TEME
    else
        rteme = rxyz; %no prec nut
    end
    
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
    
    if precnut ~= 0
        fgxyz = teme2eci(fgteme,[0;0;0],[0;0;0],ttt,106,0,'a');
    else
        fgxyz = fgteme; %no prec nut
    end
    fgtsw = xyz2tsw(fgxyz,Om,in,u);
    Tg = fgtsw(1);
    Sg = fgtsw(2);
    Wg = fgtsw(3);
else
    rteme = rxyz;
    recef = R3(gstime(jd))*rteme;    
    Tg = 0;
    Sg = 0;
    Wg = 0;
end

%nonspherics simplified (J2 Perturbation)
if Zj2 ~= 0 && Zg ==0
    Tj2 = -mu*del/r^4 * sin(2*u) * (sin(in))^2;
    Sj2 = mu*del/r^4 * (3*(sin(u))^2*(sin(in))^2-1);
    Wj2 = mu*del/r^4 * sin(u) * sin(2*in);
else
    Tj2 = 0;
    Sj2 = 0;
    Wj2 = 0;
end

%ATMOSPHERIC DRAG
if Zatm ~=0
    [ F81,F107, Kp ] = getsolar( jd, 3 );
    rho = atmgost2(recef, r-6378.14,DOY,gstime(jd),rtascsun,declsun,F107,F81,Kp); % GOST
    %rho = 2e-13;
    Tatm = -Cxa*rho*(v*1000)^2*Sm/(2*m) * 10^(-3);
else
    Tatm = 0;
end

%PROPULSION
thrust = 0.003* 10^(-3); % km! SPUTNIX
%Isp = 1470; %SPD70
%Isp = 1550; %SPD140
%Isp = 2660; %SPD100D
%Isp = 1600; %SPD100V
Isp = 1000; % SPUTNIX
fprop = thrust/m;
mdot = thrust/(Isp*g0);

Tprop = cos(PROPXSI)*fprop;
Wprop = -sign(cos(u))*sin(PROPXSI)*fprop;

%TOTAL PERTURBATIONS
T = Zj2*Tj2 + Zg*Tg + Zsunmoon*Tsunmoon + Zemp*Temp + Zprop*thrustflag*Tprop + Zatm*Tatm;
S = Zj2*Sj2 + Zg*Sg + Zsunmoon*Ssunmoon + Zemp*Semp;
W = Zj2*Wj2 + Zg*Wg + Zsunmoon*Wsunmoon + Zemp*Wemp + Zprop*thrustflag*Wprop;


%EQ!
Dy(1) = 2*r*(p/mu)^(1/2)*T;
Dy(2) = sqrt(p/mu)*(T*(1+r/p)*cos(u) + S*sin(u) + r/p*(T*l2-W*l1/tan(in)*sin(u)));
Dy(3) = sqrt(p/mu)*(T*(1+r/p)*sin(u) - S*cos(u) + r/p*(T*l1+W*l2/tan(in)*sin(u)));
Dy(4) = -r/(mu*p)^(1/2)*sin(u)/sin(in)*W;
Dy(5) = -r*cos(u)/(mu*p)^(1/2)*W;
Dy(6) = (mu*p)^(1/2)/r^2 + r/(mu*p)^(1/2)*sin(u)/tan(in)*W;
Dy(7) = -Zprop*thrustflag*mdot;

% 
% if jd >= 2457730
%     TESTR = [TESTR;[jd-jday(2016,7,1,12,0,0),p/(1-e^2),p,r-6378.14,rho, Tatm, Tg, Tsunmoon, Temp, Tprop]];
% end
%     
    


end
