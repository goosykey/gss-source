function [am, em, omm, Omm, inm, Mm, rm] = meanElementsP(a, e, om, Om, in, u)
%MEANELEMENTS Summary of this function goes here
%   Detailed explanation goes here

mu = 398600.44;

f = u-om;

T = 2*pi*sqrt(a^3/mu);

[R, V] = math.randv(a,e,in,Om,om,f);
forceModels = forceModelSet; % default set of force models
forceModels.atm = 'jr'; 
forceModels.grav= 'j2';
forceModels.sunmoon = 'none';
forceModels.emp= 'none';

[~,rvmi] = ode45(@(t,rv)earth324(t,rv,'forceModels', forceModels) ...
    , 0:T, [R;V;100]);

Ri = rvmi(:,1:3); ri = sqrt(sum(Ri.^2,2));
Vi = rvmi(:,4:6); vi = sqrt(sum(Vi.^2,2));

Ixyz = cross(Ri,Vi,2); % orbital momentum vector
I = sqrt(sum(Ixyz.^2,2));

Esp = vi.^2/2 - mu./ri;

ai = -mu./(2*Esp);
ei = sqrt(1-I.^2./(ai*mu));
pii = ai.*(1-ei.^2);


Omi = atan2(Ixyz(:,1),-Ixyz(:,2));
ini = acos(Ixyz(:,3)./I);
ui = atan2(Ri(:,3)./sin(ini), Ri(:,1).*cos(Omi)+Ri(:,2).*sin(Omi));

enull = e < 1e-06;

nui = atan2(sqrt(pii/mu).*dot(Ri,Vi,2),pii-Ri);
nui(enull) = ui(enull);

omi = ui-nui;

sinEi = (sqrt(1-ei.*ei) .* sin(nui)) ./ (1 + ei.*cos(nui));
cosEi = (ei + cos(nui)) ./ (1 + ei.*cos(nui));

Ei = atan2(sinEi, cosEi);
Mi = Ei - ei .* sin(Ei);

am = mean(ai);
em = mean(ei);
Omm = mean(Omi);
inm = mean(ini);
omm = mean(omi);
Mm = mean(Mi);
rm = mean(ri);




end