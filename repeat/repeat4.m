% repeat4.m         

% calculates the osculating semimajor axis required for
% a user-defined repeating ground track orbit

% numerical integration solution

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

global mu req j2 smu suncoef rkcoef

global omega lgrav mgrav jdate0 gst0 isun ccoef scoef

global oev tnode elong0 norbits

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% astrodynamic and utility constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Earth equatorial radius (kilometers)

req = 6378.14;

% Earth gravitational constant (kilometer^3/second^2)

mu = 398600.4415;

% Sun gravitational constant (kilometer^3/second^2)

smu = 132712440040.9446;

% Earth inertial rotation rate (radians/second)

omega = 7.2921151467e-5;

% angular conversion constants

dtr = pi / 180.0;

rtd = 180.0 / pi;

pi2 = 2.0 * pi;

% root-finding convergence criterion

rtol = 1.0e-8;

% initialize rkf78

rkcoef = 1;

% initialize sun ephemeris

suncoef = 1;

% read gravity model coefficients

[ccoef, scoef] = readgm('egm96.dat');

j2 = -ccoef(3, 1);

% begin simulation

clc; home;
   
fprintf('\n                      program repeat4\n');
   
fprintf('\n< required osculating semimajor axis - integrated solution >\n');

% request initial calendar date

[month, day, year] = getdate;

% request initial utc

[utc_hr, utc_min, utc_sec] = get_utc_time;

% request initial osculating orbital elements

oev = getoe([1;1;1;1;1;0]);

% compute julian date at 0 hours utc

day = day + utc_hr / 24.0 + utc_min / 1440.0 + utc_sec / 86400.0;

jdate0 = julian(month, day, year);

% compute greenwich sidereal time at 0 hours utc

gst0 = gast4(jdate0, 0.0, 1);

while(1)
    
   fprintf('\n\nplease input the number of integer orbits in the repeat cycle\n');

   norbits = input('? ');
   
   if (norbits > 0)
      break;
   end 
   
end

% request degree and order of gravity model

while(1)

   fprintf('\n\nplease input the degree of the Earth gravity model (zonals)\n');
   fprintf('(2 <= zonals <= 18)\n');

   lgrav = input('? ');

   if (lgrav >= 0)

      break;

   end

end

while(1)

   fprintf('\nplease input the order of the Earth gravity model (tesserals)\n');
   fprintf('(0 <= tesserals <= 18)\n');

   mgrav = input('? ');

   if (mgrav >= 0)

      break;

   end

end

while(1)

   fprintf('\nwould you like to include the point-mass gravity of the Sun (y = yes, n = no)\n');

   yn = lower(input('? ', 's'));

   if (yn == 'y' || yn == 'n')

      break;

   end

end

if (yn == 'y')

   isun = 1;

else

   isun = 0;

end

% set initial true anomaly to ascending node crossing

oev(6) = mod(-oev(4), 2.0 * pi);

% initial east longitude of the ascending node (radians)

elong0 = mod(oev(5) - gst0, pi2);

% define a "bracket" for the semimajor axis (kilometers)

x1 = oev(1) - 100.0;

x2 = oev(1) + 100.0;

% calculate osculating semimajor axis

clc; home;

fprintf('\n  working ...\n');

[xroot, froot] = brent('rpt4fun1', x1, x2, rtol);

sma = xroot;

% Keplerian period (seconds)

tkepler = 2.0 * pi * sma * sqrt(sma / mu);

% print results

clc; home;

fprintf('\n                      program repeat4\n');
   
fprintf('\n< required osculating semimajor axis - integrated solution >\n');

fprintf('\nsemimajor axis             %12.6f  kilometers \n', sma);

fprintf('\neccentricity               %12.8f  \n', oev(2));

fprintf('\ninclination                %12.6f  degrees \n', oev(3) * rtd);

fprintf('\nargument of perigee        %12.6f  degrees \n', oev(4) * rtd);

fprintf('\nraan                       %12.6f  degrees \n', oev(5) * rtd);

fprintf('\ntrue anomaly               %12.6f  degrees \n', oev(6) * rtd);

fprintf('\nKeplerian period           %12.6f  minutes \n', tkepler / 60.0);

fprintf('\naverage nodal period       %12.6f  minutes \n', tnode / 60.0);

fprintf('\nsolar days to repeat       %12.6f  days \n', (tnode * norbits) / 86400.0);

fprintf('\nnumber of orbits to repeat   %6i \n', norbits);

fprintf('\ndegree of gravity model        %2i \n', lgrav);

fprintf('\norder of gravity model         %2i \n', mgrav);

if (isun == 1)
    
    fprintf('\nsimulation includes the point-mass gravity of the Sun\n\n');
    
else
    
    fprintf('\nsimulation does not include the point-mass gravity of the Sun\n\n');
    
end



   
