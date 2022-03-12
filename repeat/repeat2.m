% repeat2.m       

% estimate time to repeat ground track

% numerically integrated solution

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

global mu req j2 smu suncoef rkcoef yi yfinal tiwrk

global omega lgrav mgrav jdate0 gst0 isun ccoef scoef

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

clc; home;

fprintf('\n                    program repeat2 \n\n');

fprintf(' < time to repeat ground track - integrated solution > \n');

% request initial calendar date

[month, day, year] = getdate;

% request initial utc

[utc_hr, utc_min, utc_sec] = get_utc_time;

% compute julian date at 0 hours utc

day = day + utc_hr / 24.0 + utc_min / 1440.0 + utc_sec / 86400.0;

jdate0 = julian(month, day, year);

% compute greenwich sidereal time at 0 hours utc

gst0 = gast4(jdate0, 0.0, 1);

% request orbital elements at ascending node

oev = getoe([1;1;1;1;1;0]);

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

while(1)

   fprintf('\n\nplease input the closure tolerance (degrees)\n');

   fprintf('(a value between 0.1 and 0.5 degrees is recommended)\n');

   dlong = input('? ');

   if (dlong > 0.0)
      break;
   end

end

clc; home;

fprintf('\n  working ...\n');

% true anomaly

oev(6) = mod(-oev(4), pi2);

% Keplerian period

tkepler = pi2 * sqrt(oev(1) * oev(1) * oev(1) / mu);

% determine initial eci state vector

[r, v] = orb2eci(mu, oev);

% load state vector working array

yold(1:3) = r;

yold(4:6) = v;

% step size guess (seconds)

h = 60.0;

% propagation step size

dtstep = 120.0;

% initialize "final" time

tf = 0.0;

% initial objective function

fxnew = 0.0;

% initialize accumulated number of orbits

norbits = 0;

% initial east longitude of the ascending node (radians)

elong0 = mod(oev(5) - gst0, pi2);

while (1)
    
    rflg = 0;
    
    % save current value of objective function
    
    fxold = fxnew;
    
    % set initial time to final time
    
    ti = tf;
    
    % increment final time
    
    tf = ti + dtstep;
    
    % save initial time as left side of bracket
    
    tisaved = ti;
    
    % integrate equations of motion
    
    neq = 6;
    
    tetol = 1.0d-10;
    
    ynew = rkf78('repeat_eqm', neq, ti, tf, h, tetol, yold);
    
    % compute current value of z component of position vector
    
    fxnew = ynew(3);
    
    % reset state vector
    
    yold = ynew;
    
    % check to see if the user-defined objective
    % function has been bracketed during this step
    
    if (fxnew * fxold < 0.0 && ynew(6) > 0.0)
        
        rflg = 1;
        
    end
    
    % if bracketed, find root of objective function
    
    if (rflg == 1)
        
        % load "working" time and state vector array
        % as values on right side of bracket
        
        tiwrk = tf;
        
        yi = ynew;
        
        % find root of objective function
        
        [troot, froot] = brent('rpt2func', tisaved, tiwrk, rtol);
                
        norbits = norbits + 1;
        
        % compute current east longitude of the ascending node (radians)
        
        a = mod(gst0 + omega * troot, pi2);
        
        b = sin(a);
        
        c = cos(a);
        
        ctmp(1) = c * yfinal(1) + b * yfinal(2);
        
        ctmp(2) = c * yfinal(2) - b * yfinal(1);
        
        elong = atan3(ctmp(2), ctmp(1));
        
        % check for convergence
        
        if (abs(elong - elong0) <= (dlong * dtr))
            
            break;
            
        end
        
    end
    
end

% print results

clc; home;

fprintf('\n              program repeat2 \n\n');

fprintf('< time to repeat ground track - integrated solution > \n');

fprintf('\nsemimajor axis             %12.6f  kilometers \n', oev(1));

fprintf('\neccentricity               %12.8f  \n', oev(2));

fprintf('\ninclination                %12.6f  degrees \n', oev(3) * rtd);

fprintf('\nargument of perigee        %12.6f  degrees \n', oev(4) * rtd);

fprintf('\nraan                       %12.6f  degrees \n', oev(5) * rtd);

fprintf('\ntrue anomaly               %12.6f  degrees \n', oev(6) * rtd);

fprintf('\nKeplerian period           %12.6f  minutes \n', tkepler / 60.0);

fprintf('\naverage nodal period       %12.6f  minutes \n', (troot / norbits) / 60.0);

fprintf('\nfinal delta-longitude      %12.6f  degrees \n', rtd * abs(elong - elong0));

fprintf('\nsolar days to repeat       %12.6f  days \n', troot / 86400.0);

fprintf('\nnumber of orbits to repeat   %6i \n', norbits);

fprintf('\ndegree of gravity model        %2i \n', lgrav);

fprintf('\norder of gravity model         %2i \n', mgrav);

if (isun == 1)
    
    fprintf('\nsimulation includes the point-mass gravity of the Sun\n\n');
    
else
    
    fprintf('\nsimulation does not include the point-mass gravity of the Sun\n\n');
    
end

