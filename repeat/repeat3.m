% repeat3.m       

% determine mean semimajor axis required
% for a repeating ground track orbit

% Wagner's algorithm

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

% astrodynamic and utility constants

j2 = 0.00108263;             % zonal gravity constant (nd)
req = 6378.14;               % Earth equatorial radius (kilometers)
mu = 398600.5;               % Earth gravitational constant (km^3/sec^2)
omega = 7.2921151467e-5;     % Earth inertial rotation rate (rad/sec)

pi2 = 2 * pi;                % 2 * pi
dtr = pi/180;                % degrees to radians
rtd = 180/pi;                % radians to degrees 

clc; home;

fprintf('\n                program repeat3 \n\n');
   
fprintf(' < repeating ground track - Wagner''s method > \n\n');

% request orbital elements

oev = getoe([0;1;1;1;0;0]);

ecc = oev(2);

inc = oev(3);

while(1)
    
   fprintf('\n\nplease input the number of orbits in the repeat cycle\n');

   norbits = input('? ');
   
   if (norbits > 0)
       
      break;
      
   end  
   
end

while(1)
    
   fprintf('\nplease input the number of nodal days in the repeat cycle\n');

   ndays = input('? ');
   
   if (ndays > 0)
       
      break;
      
   end
   
end

clc; home;

fprintf('\n  working ...\n');

c20 = -j2;

% calculate initial guess for semimajor axis

a0 = mu^(1/3) * ((norbits / ndays) * omega)^(-2/3);

aold = a0 * (1 - c20 * (req / a0)^2 * (4 * cos(inc)^2 ...
   - (norbits / ndays) * cos(inc) - 1));

while(1)
    
   slr = aold * (1 - ecc * ecc);

   tmp1 = mu^(1/3) * (norbits * omega / ndays)^(-2/3);

   tmp2 = (1 - 1.5 * c20 * (req / slr) ^ 2 * (1 - 1.5 * sin(inc)^2))^(2/3);

   tmp3 = (1 + c20 * (req / slr)^2 * (1.5 * (norbits / ndays) ...
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
c = req / slr;
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

% delta longitude per nodal period (radians)

dlong = tnode * (omega - raandot);

% length of nodal day (seconds)

tnday = pi2 / (omega - raandot);

% number of solar days to repeat

nsdays = tnode * norbits / 86400;

% print results

clc; home;

fprintf('\n                 program repeat3 \n\n');
   
fprintf('    < repeating ground track - Wagner''s method > \n\n');

fprintf('mean semimajor axis                %12.6f  kilometers \n', sma);

fprintf('mean eccentricity                  %12.6f \n', oev(2));

fprintf('mean inclination                   %12.6f  degrees \n', rtd * oev(3));

fprintf('mean argument of perigee           %12.6f  degrees \n', rtd * oev(4));

fprintf('mean raan                          %12.6f  degrees \n\n', rtd * oev(5));

fprintf('number of orbits to repeat         %12.6f \n', norbits);

fprintf('number of solar days to repeat     %12.6f \n\n', nsdays);

fprintf('Keplerian period                   %12.6f  minutes \n', tkepler / 60);

fprintf('nodal period                       %12.6f  minutes \n\n', tnode / 60);

fprintf('length of nodal day                %12.6f  minutes \n', tnday / 60);

fprintf('fundamental interval               %12.6f  degrees \n\n', dlong * rtd);

