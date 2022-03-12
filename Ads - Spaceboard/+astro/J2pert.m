function [ EPHEMERIS ] = J2pert( INITIAL_STATE, t )
%J2 perturbation propagator
%   - including symbolic math toolbox
%           by Alexander Kharlan, last rev. 10.01.2019
%   
%   
%   INITIAL STATE of the satellite (keplerian):
%       [a e om Om in u]
%       a  : semi-major axis(km)
%       e  : eccentricity
%       om : arg. of periapsis (radians)
%       Om : rt.asc. ascending node (radians)
%       in : inclination (radians)
%       u  : arg. of latitude (radians)
%
%       t  : time vector in seconds; e.g. [0, 15, 30, 45 ... 86400]

[M,N] = size(t);

if M == 1 && N ~= 1
    t = t';
end

R0 = 6378.14;
mu0 = 398600.44;
J2 = 1.7555e10;

ini = INITIAL_STATE;

a0 = ini(1);
e0 = ini(2);
om0 = ini (3);
Om0 = ini(4);
in0 = ini(5);
u0 = ini(6);

p0 = a0*(1-e0^2);
nu0 = u0-om0;
T = 2*pi*sqrt(a0^3/mu0);
n = 2*pi/T;

syms('eccan0');

if e0 == 0
    E0 = nu0;
    M0 = nu0;
else
    E0 = solve(tan(nu0)==sqrt(1-e0^2)*sin(eccan0)/(cos(eccan0)-e0),eccan0);
    M0 = E0 - e0*sin(E0);
    M0 = double(M0);
end

% E0
% M0

EPHEMERIS = zeros(numel(t),6);
EPHEMERIS(:,1) = a0;
EPHEMERIS(:,2) = e0; e = EPHEMERIS(:,2);
EPHEMERIS(:,5) = in0;



EPHEMERIS(:,3) = om0-2*pi*J2/mu0/p0^2 * 3*(5/4*sin(in0)^2-1) .*t/T;
om = EPHEMERIS(:,3);
EPHEMERIS(:,4) = Om0-2*pi*J2/mu0/p0^2 * 3/2*cos(in0) .*t/T;

% size(M0)
% size(n)
% size(t)
M = M0 + n.*t;

E = M;
eps = inf;

while max(eps) > 1e-5
    E1 = M+e.*sin(E);
    eps = abs(E-E1);
    E = E1;
end

nu = atan2( sqrt(1-e.^2).*sin(E),(cos(E)-e) );

% nu(abs(mod(nu,2*pi)-mod(E,2*pi))>pi/2) = nu(abs(mod(nu,2*pi)-mod(E,2*pi))>pi/2) + pi;
% nu = mod(nu,2*pi);

EPHEMERIS(:,6) = nu + om;
EPHEMERIS(:,6) = mod(EPHEMERIS(:,6),2*pi);




end

