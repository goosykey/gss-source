function [ EPHEMERIS ] = j2circ( INITIAL_STATE, t )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

R0 = 6378.14;
mu0 = 398600.44;
J2 = 1.7555e10;

ini = INITIAL_STATE;

Nsat = numel(ini(:,1));

a0 = ini(:,1);
e0 = ini(:,2);
om0 = ini (:,3);
Om0 = ini(:,4);
in0 = ini(:,5);
u0 = ini(:,6);

if any(e0 ~= 0)
    error('Orbit is not circular');
end

p0 = a0.*(1-e0.^2);
T = 2*pi*sqrt(a0.^3/mu0);
n = 2*pi./T;

Nt = numel(t);

EPHEMERIS = zeros(Nsat,Nt,6);
EPHEMERIS(:,:,1) = repmat(a0,[1,Nt]);
EPHEMERIS(:,:,2) = repmat(e0,[1,Nt]);
EPHEMERIS(:,:,5) = repmat(in0,[1,Nt]);

EPHEMERIS(:,:,3) = 0;

% a = size(repmat(Om0,[1,Nt]))
% b = size(repmat(p0,[1,Nt]))
% c = size()

EPHEMERIS(:,:,4) = repmat(Om0,[1,Nt]) -...
    2*pi*J2./mu0./repmat(p0,[1,Nt]).^2 .* 3/2.*repmat(cos(in0),[1,Nt]) ...
    .* repmat(t,[Nsat,1])./repmat(T,[1,Nt]);
EPHEMERIS(:,:,6) = repmat(u0,[1,Nt])+ repmat(n,[1,Nt]) .* repmat(t,[Nsat,1]);

EPHEMERIS(:,:,6) = mod(EPHEMERIS(:,:,6),2*pi);


end

