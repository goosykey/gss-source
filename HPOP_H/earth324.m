function Dy = earth324(t, y, varargin)
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
f_thrustVector = @(jd, y)[1;0;0];
f_Cxa = @(jd, y)2.2;

f_solar = @(jd, sigmaLevel)astro.getsolar(jd, sigmaLevel); % must return F107, F81, Kp

forceModels = forceModelSet; % default value of class;

precnut = 1; % always use precnut

t

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
        case 'f_solar'
            f_solar = varargin{2};
        case 'forcemodels'
            forceModels = varargin{2};
        otherwise
            %varargin{1}
            error(['Unexpected option: ' varargin{1}])
    end
    varargin(1:2) = [];
end

%% INITIALIZE

Dy = y(:);

rxyz = [y(1), y(2), y(3)]';
vxyz = [y(4), y(5), y(6)]';
m = y(7);

r = norm(rxyz);
v = norm(vxyz);

Ixyz = cross(rxyz,vxyz); % orbital momentum vector
I = norm(Ixyz);

Om = atan2(Ixyz(1),-Ixyz(2));
in = acos(Ixyz(3)/I);
u = atan2(rxyz(3)/sin(in), rxyz(1)*cos(Om)+rxyz(2)*sin(Om));

jd = jd0 + t/86400;
[month, day, ~] = astro.gdate (jd);
DOY = 30*(month-1) + round(day);


%% CONSTANTS
del = 66.07354399E03; % simplified J2 (all constants packed)
mu = 398600.4415; % Earth grav.param. km3/s2
muM = 7.36E+22 * 6.67384E-20; % Moon grav.param. km3/s2
muS = 1.98892E+30 * 6.67384E-20; % Sun grav.param. km3/s2
au1 = 149597871; % 1 Astronomic Unit in kms
R0 = 6378.14; % mean Earth radius
g0 = 9.80665 * 10^(-3); % g0 km/s2
I0 = 1366; %W/m^2 Integral solar flux at Earth orbit (Sun constant)
c = 3E+08; %m/s Speed of Light


%% SUNMOON

[rsun,rtascsun,declsun] = astro.sun ( jd ); %Vallado

if ~strcmpi(forceModels.sunmoon, 'none')
        
    r12 = astro.moon(jd)'*R0; r12abs = sqrt(sum(r12 .* r12));
    r13 = rsun*au1; r13abs = sqrt(sum(r13 .* r13));
    r2 = rxyz - r12;    r2abs = sqrt(sum(r2 .* r2));
    r3 = rxyz - r13;    r3abs = sqrt(sum(r3 .* r3));
    
    %XYZ cartesian 3body perturbations
    fsunxyz =  -muS/r3abs^3*r3 - muS/r13abs^3*r13; 
    fmoonxyz = -muM/r2abs^3*r2 - muM/r12abs^3*r12;
    if strcmp(forceModels.sunmoon,'sun')
        fsunmoon = fsunxyz;
    elseif strcmp(forceModels.sunmoon,'moon')
        fsunmoon = fmoonxyz;
    else
        fsunmoon = fsunxyz + fmoonxyz;
    end
else
    fsunmoon = [0,0,0]';   
end

%% DRAG & EMP PARAMETERS

Sm = f_Sm(jd, y);
Semp = f_Se(jd, y);
Cxa = f_Cxa(jd, y);

%% EM PRESSURE

los = astro.sight(rxyz,rsun*au1,'e');

if ~strcmpi(forceModels.emp, 'none') && strcmpi(los,'yes')
    pem = I0/c*2;
    fempxyz = -pem*Semp/m*rsun/norm(rsun')/10^3; %<<km!
else 
    fempxyz = [0;0;0];
end




%% nonspherics

if precnut ~= 0
    ttt= (jd - 2451545.0  )/ 36525.0; %Julian Cent
    rteme = astro.eci2teme(rxyz,[0;0;0],[0;0;0],ttt,106,0,'a'); %rvec J2000 to TEME
else
    rteme = rxyz; %no prec nut
end

recef = math.R3(astro.gstime(jd))*rteme;

if strcmpi(forceModels.grav, 'full') % full grav
    fgxyz = grav (recef, jd, ...
        'tides', forceModels.tides, ...
        'maxdeg', forceModels.maxdeg, 'maxord', forceModels.maxord, ...
        'rsun', rsun*au1, 'rmoon', astro.moon(jd)'*R0);
    
elseif strcmpi(forceModels.grav, 'j2') % simplified j2
    Tj2 = -mu*del/r^4 * sin(2*u) * (sin(in))^2;
    Sj2 = mu*del/r^4 * (3*(sin(u))^2*(sin(in))^2-1);
    Wj2 = mu*del/r^4 * sin(u) * sin(2*in);
    fgxyz = math.tsw2xyz([Tj2;Sj2;Wj2],Om,in,u);
else
    fgxyz = [0;0;0];
end


%% ATMOSPHERIC DRAG
if ~strcmpi(forceModels.atm, 'none')
    [ F107,F81, Kp ] = f_solar( jd, 3 );
    
    if strcmpi(forceModels.atm, 'gost')
        rho = astro.atmgost2(recef, r-6378.14,DOY,astro.gstime(jd),rtascsun,declsun,F107,F81,Kp); % GOST
    elseif strcmpi(forceModels.atm, 'j77')
        rho = astro.atmJ77(jd,rxyz,F107,F81,Kp); % Jacchia 1977
    elseif strcmpi(forceModels.atm, 'jr')
        rho = astro.atmJR(jd,rxyz,F107,F81,Kp); % Jacchia-Roberts
    else % just in case
        rho = astro.atmgost2(recef, r-6378.14,DOY,astro.gstime(jd),rtascsun,declsun,F107,F81,Kp); % GOST
    end

    Tatm = -Cxa*rho*(v*1000)^2*Sm/(2*m) * 10^(-3);
    fatm = math.tsw2xyz([Tatm; 0; 0],Om,in,u);
else
    fatm = [0;0;0];
end

%% PROPULSION

if ~isnan(thruster)
    thrust = thruster.thrust;
    Isp = thruster.Isp;
else
    thrust = 50; % was 0.3!!!
    Isp = 1550;
end


thrustflag = f_thrustflag(jd, y);
if thrustflag
    thrustVector = math.tsw2xyz(f_thrustVector(jd, y),Om,in,u);
else
    thrustVector = [0,0,0];
end

thrust = thrust * 10^(-3); % m to km
fprop = thrust/m;
mdot = -thrust/(Isp*g0);

fprop_xyz = fprop*thrustVector;


%% TOTAL PERTURBATIONS

fx = fgxyz(1) + fsunmoon(1) + fempxyz(1) + thrustflag*fprop_xyz(1) + fatm(1);
fy = fgxyz(2) + fsunmoon(2) + fempxyz(2) + thrustflag*fprop_xyz(2) + fatm(2);
fz = fgxyz(3) + fsunmoon(3) + fempxyz(3) + thrustflag*fprop_xyz(3) + fatm(3);


%% FINALLY -- EQUATIONS

Dy(1) = vxyz(1);
Dy(2) = vxyz(2);
Dy(3) = vxyz(3);

Dy(4) = -mu/r^3*rxyz(1) + fx;
Dy(5) = -mu/r^3*rxyz(2) + fy;
Dy(6) = -mu/r^3*rxyz(3) + fz;

Dy(7) = thrustflag*mdot;   


end


function fgxyz = grav (recef, jd, varargin)

%% CONST

R0 = 6378.137;
mu = 398600.4415;
muM = 7.36E+22 * 6.67384E-20; % Moon grav.param. km3/s2
muS = 1.98892E+30 * 6.67384E-20; % Sun grav.param. km3/s2
au1 = 149597871; % 1 Astronomic Unit in kms

%% VARARGIN DEFAULTS

maxdeg = 10;
maxord = 10;
precnut = true;
tides = true;
useAerospace = false;

rsun = [0,0,0];
rmoon = [0,0,0];


%% PARSE VARARGIN

while ~isempty(varargin)
    switch lower(varargin{1})
        case 'maxdeg'
            maxdeg = varargin{2};
        case 'maxord'
            maxord = varargin{2};
        case 'precnut'
            precnut = varargin{2};
        case 'tides'
            tides = varargin{2};
        case 'useaerospace' % this will omit modelling tides
            useAerospace = varargin{2};
        case 'rsun'
            rsun = varargin{2};
        case 'rmoon'
            rmoon = varargin{2};
        otherwise
            %varargin{1}
            error(['Unexpected option: ' varargin{1}])
    end
    varargin(1:2) = [];
end

if tides
    if norm(rsun) < 0.01*au1 || norm(rmoon) < 1e5
        tides = false;
        warning('Sun/moon position incorrect / not supplied for tidal modelling. Tidal modelling turned off.');
    end
end

%% IMPLEMENTATION

ttt= (jd - 2451545.0  )/ 36525.0; %Julian Cent

if tides
    rsun_ecef = math.R3(astro.gstime(jd))*rsun;
    rmoo_ecef = math.R3(astro.gstime(jd))*rmoon;
    r1sun = norm(rsun_ecef);
    r1moo = norm(rmoo_ecef);
    FI_sun = asin(rsun_ecef(3)/r1sun);
    FI_moo = asin(rmoo_ecef(3)/r1moo);
    LA_sun = atan2(rsun_ecef(2),rsun_ecef(1));
    LA_moo = atan2(rmoo_ecef(2),rmoo_ecef(1));    
end


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
            fuck = (-1)^w; %KONDON FACTOR EXCLUSION
            
            Cqw = astro.gcons('C',q,w);
            Sqw = astro.gcons('S',q,w);
            
            % corrections to C and S for tides
            % tides modelling per <PASTE REFERENCE>
            if tides
                if q == 2 || q == 3
                    K = astro.loveNumbers(q,w); % get Love numbers
                    N = astro.Nlm(q,w); % get normalization coefficient
                    leg_sun = legendre(q,sin(FI_sun)); % non-normalized Legendre with Kondon
                    leg_moon = legendre(q,sin(FI_moo));
                    legS = fuck * N * leg_sun(w+1); % normalized Legendre with no Kondon
                    legM = fuck * N * leg_moon(w+1);
                    knm_complex = K(3) + 1i * K(4); % complex Love number for Anelsatic Earth
                    delta_complex = knm_complex/(2*q+1) * ...
                        ( muS/mu * (R0/r1sun)^(q+1) * legS * exp(-1i*w*LA_sun) ...
                        + muM/mu * (R0/r1moo)^(q+1) * legM * exp(-1i*w*LA_moo));
                    dC_norm = real(delta_complex);
                    dS_norm = imag(delta_complex);
                    %fprintf('%g, %g >>> % 13.11f || % 13.11f \n', q, w, Cqw/N, Cqw/N + dC_norm);
                    Cqw = N * (Cqw/N + dC_norm);
                    Sqw = N * (Sqw/N + dS_norm);
                elseif q == 4 && (w == 0 || w == 1 || w == 2)
                    K = astro.loveNumbers(2,w);
                    N = astro.Nlm(2,w);
                    leg_sun = legendre(2,sin(FI_sun));
                    leg_moon = legendre(2,sin(FI_moo));
                    legS = fuck * N * leg_sun(w+1);
                    legM = fuck * N * leg_moon(w+1);
                    knm_complex = K(5);
                    delta_complex = knm_complex/5 * ...
                        ( muS/mu * (R0/r1sun)^3 * legS * exp(-1i*w*LA_sun) ...
                        + muM/mu * (R0/r1moo)^3 * legM * exp(-1i*w*LA_moo));
                    dC_norm = real(delta_complex);
                    dS_norm = imag(delta_complex);
                    %fprintf('%g, %g >>> % 13.11f || % 13.11f \n', q, w, Cqw/N, Cqw/N + dC_norm);
                    Cqw = N * (Cqw/N + dC_norm);
                    Sqw = N * (Sqw/N + dS_norm);                    
                end
            end
            
            % Spherical harmonic derivatives see Vallado p.550
            uro = (R0/RO)^q*(q+1)*fuck*LEG(q-1,w+1)*(Cqw*cos(w*LA)+Sqw*sin(w*LA));
            ufi = (R0/RO)^q*(fuck*(-1)*LEG(q-1,w+2)-mtP(q+1,w+1))*(Cqw*cos(w*LA)+Sqw*sin(w*LA));
            ula = (R0/RO)^q*w*fuck*LEG(q-1,w+1)*(Sqw*cos(w*LA)-Cqw*sin(w*LA));
            
            dUro = dUro - mu/RO^2*uro;
            dUfi = dUfi + mu/RO*ufi;
            dUla = dUla + mu/RO*ula;
        end
    end
    %fprintf('-----\n');
    
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