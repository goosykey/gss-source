function fx = rpt2func(x)

% rz objective function

% required by repeat2.m and repeat4.m

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global tiwrk yi yfinal

% number of differential equations

neq = 6;

% step size guess

h = 60.0;

% integrate from tiwrk to requested time = x

tetol = 1.0e-10;

yfinal = rkf78('repeat_eqm', neq, tiwrk, x, h, tetol, yi);

% compute current value of z component of position vector

fx = yfinal(3);




