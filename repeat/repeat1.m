% repeat1.m      

% time to repeat ground track

% Kozai analytic orbit propagation

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

% global constants

global j2 req mu mm pmm apdot raandot

% astrodynamic and utility constants

j2 = .00108263;              % zonal gravity constant (nd)
req = 6378.14;               % Earth equatorial radius (kilometers)
mu = 398600.5;               % Earth gravitational constant (km^3/sec^2)
omega = 7.2921151467e-5;     % Earth inertial rotation rate (rad/sec)

pi2 = 2 * pi;                % 2 * pi
dtr = pi/180;                % degrees to radians
rtd = 180/pi;                % radians to degrees 

clc; home;

fprintf('                   program repeat1 \n\n');
   
fprintf(' < time to repeat ground track - analytic solution > \n\n');

% request orbital elements

oev1 = getoe([1;1;1;1;0;0]);

while(1)

   fprintf('\n\nplease input the closure tolerance (degrees)\n');
   
   fprintf('(a value between 0.1 and 0.5 degrees is recommended)\n');

   tol = input('? ');
   
   if (tol > 0)
       
      break;
      
   end  
   
end

% set true anomaly to ascending node

oev1(6) = -oev1(4);

% calculate initial mean anomaly (radians)

a = sqrt(1 - oev1(2) * oev1(2)) * sin(oev1(6));

b = cos(oev1(6)) + oev1(2);

eanom = atan3(a, b);

oev1(6) = mod(eanom - oev1(2) * sin(eanom), 2.0 * pi);

% compute Kozai mean motion and perturbations

t = 0;

[oev2, r, v] = kozai1(1, t, oev1);

% Keplerian period (seconds)

tkepler = pi2 / mm;

% nodal period - time between nodal crossings

tnode = 2 * pi / (apdot + pmm);

% delta longitude per nodal period (radians)

dlong = tnode * (omega - raandot);

% length of nodal day (seconds)

tnday = pi2 / (omega - raandot);

% initialize

xlong = 0;

norbits = 0;

% increment nodal crossing and look for closure

while (1)
    
  % increment current longitude
     
  xlong = xlong + dlong;
  
  % mod longitude

  if (xlong >= pi2)
      
     xlong = xlong - pi2;
     
  end
     
  % increment number of orbits

  norbits = norbits + 1;
     
  % check for ground track closure

  if (abs(xlong - pi2) <= dtr * tol || abs(xlong) <= dtr * tol)
     break;
  end
  
end

% number of days to repeat ground track

ndays = tnode * norbits / 86400;

% actual closure delta longitude (radians)

delta = abs(xlong - pi2);
  
if (abs(xlong) < delta) 
   delta = abs(xlong);
end

% print results

clc; home;

fprintf('\n\n                    program repeat1 \n\n');
   
fprintf('  < time to repeat ground track - analytic solution > \n\n');

fprintf('mean semimajor axis                %12.6f  kilometers \n', oev1(1));

fprintf('mean eccentricity (nd)             %12.6f \n', oev1(2));

fprintf('mean inclination                   %12.6f  degrees \n', rtd * oev1(3));

fprintf('mean argument of perigee           %12.6f  degrees \n', rtd * oev1(4));

fprintf('mean raan                          %12.6f  degrees \n\n', rtd * oev1(5));

fprintf('number of orbits to repeat         %12.6f \n', norbits);

fprintf('number of solar days to repeat     %12.6f \n\n', ndays);

fprintf('Keplerian period                   %12.6f  minutes \n', tkepler / 60);

fprintf('nodal period                       %12.6f  minutes \n\n', tnode / 60);

fprintf('length of nodal day                %12.6f  minutes \n', tnday / 60);

fprintf('fundamental interval               %12.6f  degrees \n\n', dlong * rtd);

fprintf('closure tolerance                  %12.6f  degrees \n', tol);

fprintf('actual closure                     %12.6f  degrees \n\n', delta * rtd);

