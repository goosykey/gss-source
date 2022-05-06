function [Sm] = Sm_idealPanels(jd, y)
%SM_IDEALPANELS CSArea for ideal solar panels
%   This function returns cross-section area for atmospheric drag (along
%   velocity vector) for a satellite with solar panels that always point at
%   the Sun (by some ideal mechanical structure).
% 
%   INPUT:
%       jd : epoch in Julian Day (non-modified)
%       y  : phase vector in earth353 format:
%          [p
%           lambda1
%           lambda2
%           Omega
%           inc
%           u]
%
%   OUTPUT:
%       Sm : total cross-section area, in sq.m


%% INITIALIZE

S_body = 1; % m^2, area of the satellite body (not rotating)
S_batt = 1; % m^2, this area will always point at the Sun

%% CONSTANTS

au1 = 149597871; % 1 Astronomic Unit in kms


%% IMPLEMENTATION

Y = math.unorbitize(y);

a = Y(1);
e = Y(2);
om = y(3);
Om = y(4);
in = Y(5);
u = Y(6);

nu = u - om;

[~,V] = randv(a,e,in,Om,om,nu);

[rsun,~,~] = astro.sun ( jd ); %Vallado

rsun = rsun * au1;

N_dir = rsun/norm(rsun);
V_dir = V/norm(V);

Sm = S_body + S_batt * dot(N_dir, V_dir);


end