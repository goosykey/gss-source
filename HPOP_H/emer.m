clear; clc;

%% CONST

R0 = 6378.14;

deg = pi/180;

%% INIT

orb = [R0+600 0 0 0 astro.sunsyn(600)*deg 0 ];

[R0,V0] = math.randv(orb(1),orb(2),orb(5),orb(4),orb(3),orb(6)-orb(3)); % Convert to cartesian

jd0 = astro.jday(2020, 10, 23, 12, 0, 0); % Julian day

tspan = 0:60:2000;
options = odeset('MaxStep', 30,'RelTol', 1e-05);



%% EARTH324 SOLUTION


PARAMETERS = {'jd0', jd0};

[~, yi] = ode45( @(t,rv)earth324(t,rv,PARAMETERS{:}),...
    tspan, [R0;V0;200],options);

vis_earthdraw(jd0);
plot3(yi(:,1),yi(:,2),yi(:,3),'r','linewidth',2);

[~, ~, ~, Om, in, u] = math.cart2kep(yi(end,1),yi(end,2),yi(end,3),yi(end,4),yi(end,5),yi(end,6));

v_tsw = math.xyz2tsw(yi(end,4:6)',Om, in, u) + [0;0;-4];

R0 = yi(end,1:3); R0 = R0';
V0 = math.tsw2xyz(v_tsw,Om,in,u);

[~, yi] = ode45( @(t,rv)earth324(t,rv,PARAMETERS{:}),...
    tspan/2, [R0;V0;200],options);

plot3(yi(:,1),yi(:,2),yi(:,3),'c','linewidth',2);

R0 = yi(end,1:3); R0 = R0';
V0 = yi(end,4:6); V0 = V0';
PARAMETERS = {'jd0', jd0, 'f_thrustflag', @(t,y)true} ;

[~, yi] = ode45( @(t,rv)earth324(t,rv,PARAMETERS{:}),...
    tspan, [R0;V0;200],options);

plot3(yi(:,1),yi(:,2),yi(:,3),'y','linewidth',2)
 % 
HPOP/HPOP_DEMO.m
 
clear; clc;

%% ADD PATH

mfn = mfilename('fullpath');

if isunix
    slash = '/';
else
    slash = '\';
end


while mfn(end) ~= slash
    mfn(end) = [];
end

addpath(genpath([mfn slash '..' slash '..' slash 'MATLAB' slash]))

%% CONST

R0 = 6378.14;

deg = pi/180;

%% INIT

orb = [R0+600 0 0 0 astro.sunsyn(600)*deg 0 % 600 km Sun Sync
    R0+390 0 0 0 51*deg 0 % ISS orbit
    R0+35786 0 0 0 0 0 % GEO
    R0+1200 0.1 0 0 0 0 % MEO elliptical
    26553.4 0.740969 270*deg 123.375*deg 63.4*deg 0 % Molniya
    R0+220 0 0 0 90*deg 0   % Very low orbit
    R0+220 0 0 0 90*deg 0]; % Same as above but will use thruster

[R0,V0] = math.randv(orb(:,1),orb(:,2),orb(:,5),orb(:,4),orb(:,3),orb(:,6)-orb(:,3)); % Convert to cartesian

jd0 = astro.jday(2020, 10, 23, 12, 0, 0); % Julian day

tspan = 0:60:86400;
options = odeset('MaxStep', 30,'RelTol', 1e-05);



%% EARTH324 SOLUTION

RESULTS = zeros(7, numel(tspan), 6); % pre-allocate memory for results

forceModels0 = forceModelSet; % default set of force models
forceModels0.atm = 'none'; % this will omit atmospheric modelling

forceModels1 = forceModelSet;
forceModels1.atm = 'j77'; % use Jacchia 1977 atmosphere

forceModels2 = forceModelSet;
forceModels2.grav = 'j2'; % use simplified j2 gravity models

forceModels3 = forceModelSet('grav','none', 'maxdeg', 0, 'maxord', 0, ...
    'atm', 'none','tides', 'none', ...
    'sunmoon', 'none', 'emp', 'none'); % no any perturbing forces


PARAMETERS = { {'jd0', jd0} % first satellite will go with built-in default force models and parameter functions
    {'jd0', jd0, 'forceModels', forceModels1} % this satellite will use Jacchia atmosphere
    {'jd0', jd0, 'forceModels', forceModels0} % this GEO satellite will use explicity-passed default force model set
    {'jd0', jd0, 'forceModels', forceModels3} % this MEO satellite will have an unperturbed orbit
    {'jd0', jd0, 'forceModels', forceModels2} % this Molniya satellite will use simplified j2 gravity
    {'jd0', jd0}
    {'jd0', jd0, 'f_thrustflag', @(t,y)true} }; % this satellite will have its thruster on at all times

for q = 1:7
    fprintf('Modelling satellite # %g... \n', q);
    [~, yi_here] = ode45( @(t,rv)earth324(t,rv,PARAMETERS{q}{:}),...
        tspan, [R0(:,q);V0(:,q);200],options);
    RESULTS(q,:,:) = permute(yi_here(:,1:6),[3,1,2]);
end


vis_movie(tspan, RESULTS);




