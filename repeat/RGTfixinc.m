% RGTfixinc.m       July 10, 2015

% determine mean semimajor axis required
% for a repeating ground track orbit - GIVEN INCLINATION
% This is a single step in RGTSSO.m loop - to derive SUN SYN

% Wagner's algorithm - ITERATION LOOP

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [a, in] = RGTfixinc(norbits, ndays, incdeg)

% astrodynamic and utility constants

j2 = 0.0010826267;             % zonal gravity constant (nd)
R0 = 6378.14;               % Earth equatorial radius (kilometers)
mu = 398600.44;               % Earth gravitational constant (km^3/sec^2)
omega = 7.2921151467e-5;     % Earth inertial rotation rate (rad/sec)

pi2 = 2 * pi;                % 2 * pi
deg = pi/180;

%clc; home;

ecc = 0;
inc = incdeg*deg;

c20 = -j2;

% calculate initial guess for semimajor axis

a0 = mu^(1/3) * ((norbits / ndays) * omega)^(-2/3);

aold = a0 * (1 - c20 * (R0 / a0)^2 * (4 * cos(inc)^2 ...
   - (norbits / ndays) * cos(inc) - 1));

while(1)
    
   slr = aold * (1 - ecc * ecc);

   tmp1 = mu^(1/3) * (norbits * omega / ndays)^(-2/3);

   tmp2 = (1 - 1.5 * c20 * (R0 / slr) ^ 2 * (1 - 1.5 * sin(inc)^2))^(2/3);

   tmp3 = (1 + c20 * (R0 / slr)^2 * (1.5 * (norbits / ndays) ...
      * cos(inc) - 0.75 * (5 * cos(inc)^2 - 1)))^(2/3);

   anew = tmp1 * tmp2 * tmp3;

   if (abs(anew - aold) < 0.000001)
      break;
   else
      aold = anew;
   end
end

sma = anew;

% Keplerian period (seconds)

tkepler = pi2 * sma * sqrt(sma / mu);
  
% keplerian mean motion (rad/sec)

mm = sqrt(mu / (sma * sma * sma));
      
% orbital semiparameter (kilometers)

slr = sma * (1 - ecc * ecc);

b = sqrt(1 - ecc * ecc);
c = R0 / slr;
d = c * c;
e = sin(inc) * sin(inc);

% perturbed mean motion (rad/sec)

pmm = mm * (1 + 1.5 * j2 * d * b * (1 - 1.5 * e));

% argument of perigee perturbation (rad/sec)

apdot = 1.5 * j2 * pmm * d * (2 - 2.5 * e);

% raan perturbation (rad/sec)

raandot = -1.5 * j2 * pmm * d * cos(inc);

% nodal period - time between nodal crossings

tnode = 2 * pi / (apdot + pmm);

% vallado: simplified (1) and full (2)
% tnode1 = tkepler * (1 - 3*j2/2 * (R0/sma)^2 * (3-4*sin(inc)^2) );
% tnode2 = tkepler / (1 + 3*j2/4 * (R0/slr)^2 * ...
%     (sqrt(1-ecc^2)*(2-3*sin(inc)^2) + (4-5*sin(inc)^2)));

% delta longitude per nodal period (radians)

dlong = tnode * (omega - raandot);

% length of nodal day (seconds)

tnday = pi2 / (omega - raandot);

% number of solar days to repeat

nsdays = tnode * norbits / 86400;

% print results

%clc; home;

fprintf('mean semimajor axis                %12.6f  kilometers \n', sma);

fprintf('mean inclination                   %12.6f  degrees \n', 1/deg * inc);

fprintf('number of orbits to repeat         %12.6f \n', norbits);

fprintf('number of solar days to repeat     %12.6f \n\n', nsdays);

fprintf('Keplerian period                   %12.6f  minutes \n', tkepler / 60);

fprintf('nodal period                       %12.6f  minutes \n\n', tnode / 60);

fprintf('length of nodal day                %12.6f  minutes \n', tnday / 60);

fprintf('fundamental interval               %12.6f  degrees \n', dlong * 1/deg);

fprintf('new inclination                    %12.6f  degrees \n\n', astro.sunsyn(sma-R0));

a = sma;
in = astro.sunsyn(a-R0);

end

