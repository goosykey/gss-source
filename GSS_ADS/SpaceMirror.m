clear all;

Rspotreq = 25e3; % Required spot diameter

alphaMaylar = 32/60; % ???? ????????? ?????????
ro = 0.92; % mirror reflectance Maylar or Kapton coated with aluminum [1,2]
solarflux = 1360;  % Solar flux intensity, W/m2
h = 550e3; % altitude of the circle orbit, m

Rspotreal = h*tand(alphaMaylar/2);

alphareq = radtodeg(2*atan(Rspotreq/h));

mreq = -8;
Irequired = 10.^(-2*mreq/5)*2.54e-6/683; % Required Solar flux [W/m2]

tetta = linspace(30,80,100);
gamma = linspace(20,120,100)';

tau = 0.1283 + 0.7559*exp(-0.3878*secd(90-tetta)); % atmospheric transmissivity [3]
Srefl = (Irequired*pi*(h*tand(alphareq/2))^2)./(solarflux*ro*tau.*abs(cosd(gamma/2)).*abs(sind(tetta)));
Drefl = 2*sqrt(Srefl/pi);

surf(tetta, gamma, Drefl);

surf(tetta, gamma, Srefl);
Dmax = max(max(Drefl));
Smax = max(max(Srefl));
Sreflfinal = pi*(ceil(Dmax)/2)^2;

mPossible = GetMagnitude('sm', Sreflfinal, 'tetta', tetta, 'gamma', gamma, 'alpha', alphareq);

mPossiblemin = max(max(mPossible));
surf (tetta, gamma, mPossible);
