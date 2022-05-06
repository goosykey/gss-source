function [Sm] = Sm_1axisPanels(jd, y)
%Sm_1axisPanels CSArea for 1-axis solar panels
%   This function returns cross-section area for atmospheric drag (along
%   velocity vector) for a satellite with solar panels that rotate only
%   around the orbit normal to maximize the harvested solar energy.
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
S_batt = 1; % m^2, this area will rotate around W axis to maximize power

%% CONSTANTS

au1 = 149597871; % 1 Astronomic Unit in kms


%% IMPLEMENTATION

Y = y;

Y(2) = sqrt(y(2).^2 + y(3).^2);
Y(3) = atan2 (y(2),y(3));

Y(1) = y(1) ./ (1-Y(2).^2);

a = Y(1);
e = Y(2);
om = y(3);
Om = y(4);
in = Y(5);
u = Y(6);

nu = u - om;

[R,V] = math.randv(a,e,in,Om,om,nu);

[rsun,~,~] = astro.sun ( jd ); %Vallado

rsun = rsun * au1;
sun_dir = rsun/(norm(rsun));

Wminus = cross(R,V);
Wminus = Wminus/norm(Wminus);

sun_dir_tsw= math.xyz2tsw(sun_dir,Om,in,u);

n = [sun_dir_tsw(1), sun_dir_tsw(2), 0];
n = n/norm(n);

V_dir = V/norm(V);
V_dir_tsw = math.xyz2tsw(V_dir,Om,in,u);

Sm = S_body + S_batt * dot(n, V_dir_tsw);


end