function fx = rpt4fun1 (x)

% delta-fundamental interval objective function

% required by repeat4.m

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global mu omega gst0 oev tiwrk yi yfinal tnode

global norbits elong0

% extract current value of semimajor axis (kilometers)

oev(1) = x;

% determine initial eci state vector

[r, v] = orb2eci(mu, oev);

% save state vector as "old" values

yold(1:3) = r;

yold(4:6) = v;

% event search step size guess (seconds)

h = 120.0;

% orbit propagation step size (seconds)

dtstep = 60.0;

% initialize 'final' time

tf = 0.0;

% set initial value of objective function

% since we are starting at the ascending node
% the z component of the position vector = 0

fxnew = 0.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% numerically integrate and find time of nodal crossing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% truncation error tolerance

tetol = 1.0e-10;

% number of differential equations

neq = 6;

norb = 0;

while(1)
    
   % initialize 'root found ' indicator

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

   ynew = rkf78('repeat_eqm', neq, ti, tf, h, tetol, yold);

   % save current value of objective function
   % as the z component of the position vector

   fxnew = ynew(3);
   
   yold = ynew;
   
   % check to see if a nodal crossing
   % has been bracketed during this step

   if (fxnew * fxold < 0.0 && ynew(6) > 0.0)
       
      rflg = 1;
      
   end   

   % if bracketed and this is an ascending node,
   % find time of ascending node crossing

   if (rflg == 1)

      % load 'working' time and state vector array
      % as values on right side of bracket

      tiwrk = tf;

      yi = ynew;
      
      % root-finding convergence criterion
      
      rtol = 1.0e-8;

      % find time of ascending node crossing

      [troot, froot] = brent('rpt2func', tisaved, tiwrk, rtol);
            
      norb = norb + 1;
      
   end
   
   if (norb == norbits)
       
       break;
       
   end
   
end

tnode = (troot / norbits);

% compute current east longitude of the ascending node (radians)

a = mod(gst0 + omega * troot, 2.0 * pi);

b = sin(a);

c = cos(a);
      
c1 = c * yfinal(1) + b * yfinal(2);

c2 = c * yfinal(2) - b * yfinal(1);
      
elong = atan3(c2, c1);

% compute objective function as difference between current
% earth longitude and initial earth longitude of ascending node

fx = elong0 - elong;


