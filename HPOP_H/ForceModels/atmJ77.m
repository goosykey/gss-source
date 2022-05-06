function [rho] = atmJ77(jd, Rxyz, F107, F81, kp)
%ATMJ77 Summary of this function goes here
%   Detailed explanation goes here

%% CONST

R0 = 6378.137;
deg = pi/180;

f = 1/298.257223563; % ellipsoid flattening
e2 = 2*f-f^2; % ellipsoid first eccentricity squared

Av = 6.02257e23; % Avogadro





%% LOCAL COORDS

rI = Rxyz(1);
rJ = Rxyz(2);
rK = Rxyz(3);

%% SUN

[rsun, ~, declsun] = astro.sun(jd);

rx = rsun(1);
ry = rsun(2);


%% EVALUATING TEMPERATURE

% Nighttime global exospheric temperature
T_c = 379 + 3.24 * F81 + 1.3 * (F107-F81);

% coefs from Vallado p. 1001
LHA_sun = (rx*rJ - ry*rI) / abs(rx*rJ - ry*rI) ...
    * acos( (rx*rI + ry*rJ)/sqrt(rx^2+ry^2)/sqrt(rI^2+rJ^2) );
fi_gd = atan(1/(1-f)^2 * rK / sqrt(rI^2 + rJ^2));
eta = abs(fi_gd - declsun) / 2;
theta = abs(fi_gd + declsun) / 2;
tau = LHA_sun - 37*deg + 6*deg * sin(LHA_sun + 43*deg);

%% ELLIPSOIDAL PARAMETERS

% auxiliary quantities (Vallado p.138 3-7)
C0 = R0/sqrt(1-e2*sin(fi_gd)^2);
S0 = R0*(1-e2)/sqrt(1-e2*sin(fi_gd)^2);
b = (C0*cos(fi_gd)^2 + S0*sin(fi_gd)^2);
c = C0*C0*cos(fi_gd)^2 + S0*S0*sin(fi_gd)^2 - norm(Rxyz)^2;
h_ellp = -b + sqrt(b^2 - c);

%% EVALUATING TEMPERATURE (continues)

% uncorrected Exospheric temperature
T_unc = T_c * (1 + 0.3 * (sin(theta)^2.2 + (cos(eta)^2.2-sin(theta)^2.2)*cos(tau/2)^3));

deg1 = 1;

if h_ellp > 200
    dT_corr = 28*deg1 * kp + 0.03*exp(kp);
else
    dT_corr = 14*deg1 * kp + 0.02*exp(kp);
end

% corrected exospheric temperature
T_corr = T_unc + dT_corr;

%% DENSITY FROM JACCHIA 77

[~, ~, ~, ~, ~, ~, ~, CM, WM] = astro.j77sri(h_ellp, T_corr, h_ellp);

rho = CM(end)/Av .* WM(end);

rho = rho*1e3;


end