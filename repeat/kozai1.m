function [oev2, r, v] = kozai1 (iniz, t, oev1)

% analytic orbit propagation

% Kozai's method - ECI version

% input

%  iniz = initialization flag (1 = initialize, 0 = bypass)
%  t    = elapsed simulation time (seconds)

%  initial orbital elements

%  oev1(1) = semimajor axis (kilometers)
%  oev1(2) = orbital eccentricity (non-dimensional)
%            (0 <= eccentricity < 1)
%  oev1(3) = orbital inclination (radians)
%            (0 <= oev1(3) <= pi)
%  oev1(4) = argument of perigee (radians)
%            (0 <= oev1(4) <= 2 pi)
%  oev1(5) = right ascension of ascending node (radians)
%            (0 <= oev1(5) <= 2 pi)
%  oev1(6) = mean anomaly (radians)
%            (0 <= oev1(6) <= 2 pi)

% output

%  final orbital elements and state vector at time = t

%  oev2(1) = semimajor axis (kilometers)
%  oev2(2) = orbital eccentricity (non-dimensional)
%            (0 <= eccentricity < 1)
%  oev2(3) = orbital inclination (radians)
%            (0 <= oev2(3) <= pi)
%  oev2(4) = updated argument of perigee (radians)
%            (0 <= oev2(4) <= 2 pi)
%  oev2(5) = updated right ascension of ascending node (radians)
%            (0 <= oev2(5) <= 2 pi)
%  oev2(6) = updated mean anomaly (radians)
%            (0 <= oev2(6) <= 2 pi)
%  oev2(7) = updated true anomaly (radians)
%            (0 <= oev2(6) <= 2 pi)
%  r       = eci position vector (kilometers)
%  v       = eci velocity vector (kilometers/second)

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global j2 req mu mm pmm apdot raandot sinc cinc slr

% set final orbital elements equal to initial

oev2 = oev1;

if (iniz == 1)
    sinc = sin(oev1(3));
    cinc = cos(oev1(3));

    % keplerian mean motion

    mm = sqrt(mu / (oev1(1) * oev1(1) * oev1(1)));

    % orbital semiparameter

    slr = oev1(1) * (1 - oev1(2) * oev1(2));

    b = sqrt(1 - oev1(2) * oev1(2));
    c = req / slr;
    d = c * c;
    e = sinc * sinc;

    % perturbed mean motion

    pmm = mm * (1 + 1.5 * j2 * d * b * (1 - 1.5 * e));

    % argument of perigee perturbation

    apdot = 1.5 * j2 * pmm * d * (2 - 2.5 * e);

    % raan perturbation

    raandot = -1.5 * j2 * pmm * d * cinc;
end

% set initial values of argument of perigee,
% raan and mean anomaly

argperi = oev1(4);
raani = oev1(5);
manomi = oev1(6);

% propagate orbit

argper = mod(argperi + apdot * t, 2.0 * pi);
raan = mod(raani + raandot * t, 2.0 * pi);
manom = mod(manomi + pmm * t, 2.0 * pi);

oev2(4) = argper;
oev2(5) = raan;
oev2(6) = manom;

% solve Kepler's equation using Danby's method

ecc = oev1(2);

[ea, tanom] = kepler1(manom, ecc);

oev2(7) = tanom;

% position magnitude

rmag = slr / (1 + ecc * cos(tanom));

sargper = sin(argper);
cargper = cos(argper);

% argument of latitude

arglat = mod(argper + tanom, 2.0 * pi);

sarglat = sin(arglat);
carglat = cos(arglat);

sraan = sin(raan);
craan = cos(raan);

% eci position vector

r(1) = rmag * (craan * carglat - cinc * sarglat * sraan);
r(2) = rmag * (sraan * carglat + cinc * sarglat * craan);
r(3) = rmag * sinc * sarglat;

% eci velocity vector

c1 = sqrt(mu / slr);
c2 = ecc * cargper + carglat;
c3 = ecc * sargper + sarglat;

v(1) = -c1 * (craan * c3 + sraan * cinc * c2);
v(2) = -c1 * (sraan * c3 - craan * cinc * c2);
v(3) = c1 * c2 * sinc;


