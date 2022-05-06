function [am, em, omm, Omm, inm, Mm, rm] = meanElements(a, e, om, Om, in, u)
%MEANELEMENTS Summary of this function goes here
%   Detailed explanation goes here

f = u-om;

sinE = (sqrt(1-e.*e) .* sin(f)) ./ (1 + e.*cos(f));
cosE = (e + cos(f)) ./ (1 + e.*cos(f));

E = atan2(sinE, cosE);
M = E - e .* sin(E);

p = a*(1-e^2);
r = p / (1 + e*cos(f));

delta = ones(1,6);

fprintf('Iterating...');

while any(abs(delta) > eps)
    fprintf('|');
    [am, em, omm, Omm, inm, Mm, rm] = oneStep(a, e, om, Om, in, M, r);
    delta = [am-a, em-e, omm-om, Omm-Om, inm-in, Mm-M];
    %[am, em, omm, Omm, inm, Mm, rm]
    %delta
    %em
    a = am; e = em; om = omm; Om = Omm; in = inm; M = Mm; r = rm;
end

fprintf('\n');

end

function [am, em, omm, Omm, inm, Mm, rm] = oneStep(a, e, om, Om, in, M, r)

J2 = 1082.645e-06;
% J3 = -2.546e-06;
% J4 = -1.649e-06;

R = 6378.137;
%mu = 398600.4418;

E = M;
delta = inf;

while delta > eps
    E1 = M + e.*sin(E);
    delta = abs(E-E1);
    E = E1;
end

f = atan2( sqrt(1-e.^2).*sin(E),(cos(E)-e) );

u = om + f;

p = a*(1-e^2);


%% SHORTENING COMMON EXPRESSIONS

e2 = e^2;

A = 1.5*sin(in)^2;
B = (1-e2)^1.5;
C = J2*(R/p)^2;

if e == 0
    D = 0; % otherwise it will produce NaN
else
    D = 1/e * (1+1.5*e2-B);
end

cosf = cos(f);
sinf = sin(f);
cos2f = cos(2*f);
sin2f = sin(2*f);
cos3f = cos(3*f);
sin3f = sin(3*f);
cos2u = cos(2*u);
sin2u = sin(2*u);
cos2w = cos(2*om);
sin2w = sin(2*om);
sini2 = sin(in)^2;
sin2i = sin(2*in);
cosi = cos(in);

wf21 = 2*om + f;
wf2m1 = 2*om - f;
wf23 = 2*om + 3*f;
wf24 = 2*om + 4*f;
wf25 = 2*om + 5*f;

coswf21 = cos(wf21);
sinwf21 = sin(wf21);
coswf2m1 = cos(wf2m1);
sinwf2m1 = sin(wf2m1);
coswf23 = cos(wf23);
sinwf23 = sin(wf23);
coswf24 = cos(wf24);
sinwf24 = sin(wf24);
coswf25 = cos(wf25);
sinwf25 = sin(wf25);


%% 1ST ORDER SHORT-PERIODIC VARIATIONS

a_sp = J2*(R^2/a) * ((a/r)^3 * (1-A + A*cos2u) -(1-A)/B);

e_sp = 0.5*C * (1-A) * (D + 3*(1+e2/4)*cosf + 1.5*e*cos2f + e2/4*cos3f) ...
    +3/8*C*sini2 * ( (1+11/4*e2)*coswf21 + e2/4*coswf2m1 + 5*e*cos2u ...
    +1/3*(7+17/4*e2)*coswf23 + 1.5*e*coswf24 + e2/4*coswf25 + 1.5*e*cos2w );

i_sp = 3/8*C*sin2i* ( e*coswf21+cos2u+e/3*coswf23 );

om_sp = 3/4*C*(4-5*sini2)*(f-M+e*sinf) + 1.5*C*(1-A) ...
    * (1/e*(1-0.25*e2)*sinf + 0.5*sin2f+1/12*e*sin3f) ...
    -1.5*C*(1/e*(0.25*sini2+e2/2*(1-15/8*sini2))*sinwf21 + e/16*sini2*sinwf2m1 ...
    +0.5*(1-2.5*sini2)*sin2u -1/e*(7/12*sini2-e2/6*(1-19/8*sini2))*sinwf23 ...
    -3/8*sini2*sinwf24-1/16*e*sini2*sinwf25) -9/16*C*sini2*sin2w;

if isnan(om_sp)
    om_sp = 0;
end

Om_sp = -1.5*C*cosi*(f - M + e*sinf - e/2*sinwf21 - 0.5*sin2u - e/6*sinwf23);

M_sp = -1.5*C*sqrt((1-e2)/e) * ( (1-A) * ((1-0.25*e2)*sinf + e/2*sin2f + e2/12*sin3f) ...
    +0.5*sini2*( -0.5*(1+1.25*e2)*sinwf21 ...
    -e2/8*sinwf2m1 + 7/6*(1-e2/28)*sinwf23 + 0.75*e*sinwf24 + e2/8*sinwf25 )) ...
    +9/16*C*sqrt(1-e2)*sini2*sin2w;

if isinf(M_sp)
    M_sp = 0;
end

r_sp = -0.5*C*p*(1-A)*(1+e/(1+sqrt(1-e2))*cosf + 2/sqrt(1-e2)*r/a ) ...
    +0.25*C*p*sini2*cos2u;

%% FINALLY

am = a - a_sp;
em = e - e_sp;
omm = om - om_sp; %omm = mod(omm,2*pi);
Omm = Om - Om_sp;
inm = in - i_sp;
Mm = M - M_sp;
rm = r - r_sp;

%um = mod(omm + Mm, 2*pi)

end