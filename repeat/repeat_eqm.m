function ydot = repeat_eqm (t, y)

% first order equations of orbital motion

% required by repeat2.m and repeat4.m

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global jdate0 isun smu

ydot(1) = y(4);
ydot(2) = y(5);
ydot(3) = y(6);

% geopotential perturbations

agrav = gravity(t, y);

if (isun == 1)
    
   % solar perturbations

   jdate = jdate0 + t / 86400.0;

   [rasc, decl, rsun] = sun2 (jdate);

   rsunm = norm(rsun);

   % compute heliocentric position vector of the spacecraft

   for i = 1:1:3
       
       rs2sc(i) = y(i) - rsun(i);
       
   end

   rs2scm = norm(rs2sc);

   rrs2sc = smu / rs2scm ^ 3;

   rrsun = smu / rsunm ^ 3;

   for i = 1:1:3
       
       asun(i) = -rs2sc(i) * rrs2sc - rsun(i) * rrsun;
       
   end
   
else
    
   for i = 1:1:3
       
       asun(i) = 0.0;
       
   end
   
end

% total acceleration vector

for i = 1:1:3
    
    ydot(i + 3) = agrav(i) + asun(i);
    
end

