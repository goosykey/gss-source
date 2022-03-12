%% READ FILES && INITIALIZE DATA

tic;
cdb = readcities('cities1000.txt', 100);        % read city database
uidb = readui('UItable1.csv');                  % read usage intensity database
[ map, mapdata ] = readmap( 'glp15ag.asc' );    % read population density map
succ = false;

gd0 = [2017 03 20 0 0 0]; % SET UP STARTTIME (UTC)

vtau = 100; % specify voice session duration
dtau = 60;  % specify dataN session duration
stau = 400; % specify dataS session duration

Vminday = 0.001;
DMBmo = 60000;
SGBmo = 4;

br_m = 512e06;
br_f = 1e06;

jd = jday(gd0(1),gd0(2),gd0(3),gd0(4),gd0(5),gd0(6)); % find starttime Julian Date
deg = pi/180;

a = 6978.14; e = 0; om = 0; in = 97.792*deg; % specify initial orbit parameters
Om = 0*deg; u = 120*deg;

dOm = 11.125*deg; du = 8*deg; ddu = 12*deg;  % specify constellation parameters
netdata = {dOm, du, ddu, vtau, dtau, stau, Vminday, DMBmo, SGBmo, br_m, br_f};
Tb = toc;

% Preliminary stage finished

%% MASK CALCULATION
% REPEAT UNTIL SUCCESS : in case of unexpected 'parfor' loop bugs

tic;
while ~succ
    if exist('BIGASSMASK', 'var') && ~exist('OLD','var')
        mailbot('alexander.kharlan@yaliny.com','SKIP! The mask already exists. Proceeding to the simulation.','MATLAB MAIL BOT'); % send report
        break
    end
    
    try
        [ BIGASSMASK ] = intmasknew( map, mapdata, cdb, uidb );
        %[ BIGASSMASK ] = intmaskone( map, mapdata, cdb, uidb ); % WITHOUT parallel calculations
        succ = true;
        mailbot('alexander.kharlan@yaliny.com','The mask calculation has been successful!','MATLAB MAIL BOT');
    catch ER
        mailbot('alexander.kharlan@yaliny.com','The mask calculation has failed.','MATLAB MAIL BOT'); % send failure report
        delete(gcp('nocreate')); % shutdown parallel pool in case of failure
        rethrow(ER);
    end
end

Tm = toc;

% Mask calculation finished

[ Vmask, Dmask, Smask, Rmask, Omask, Pmask, SPmask, Nmask, Zmask] = decompose( BIGASSMASK );


%% OCEAN PEOPLE MASK

tic;
if ~exist('map_oc', 'var')
    [ map_oc ] = oceanpeople( mapdata, 5e06, 4e06, 200, 0.5, Rmask, Omask ); % Generate subscribers in ocean
end

To = toc;

%% GENERATE EPHEMERIS
% use J2 perturbation (secular only)

tic;

fof = J2pert([a e om Om in u],86400*3,20); % time period and tick are specified HERE

%% SIMULATION

try
    [ Nreal, Ntot, Ndes, Nass, POWER ] = CYCLO_SIM( fof, netdata, jd, map+map_oc, mapdata, BIGASSMASK  ); % get P_out from session distributions
    [ POWERS, SUNANGLES, SUBSAT, Rxyz ] = CYCLO_FIN( fof, jd, POWER, Nreal, map, mapdata );          % get everything else
    mailbot('alexander.kharlan@yaliny.com','The simulation has been successful! \n','MATLAB MAIL BOT'); % report simulation success
catch ER
    mailbot('alexander.kharlan@yaliny.com','The simulation has failed. \n','MATLAB MAIL BOT'); % report simulation failure
    rethrow(ER);
end

%% PLOT AND FORM RESULT TABLE FOR EXPORT

CYCLO_PLOT( gd0, fof, Ntot, Ndes, Nass, Nreal, POWERS, SUNANGLES )
t = fof(:,1);

MATR = [fof, Rxyz', SUBSAT, Ntot, Ndes, Nass, Nreal, Ndes-Nreal, POWERS, SUNANGLES];

Ts = toc;

%% E-MAIL REPORT ANS EXPORT RESULTS

save('A:\MATLABAUTO.mat');
csvwrite('RESULT.csv',MATR);

s = whos('BIGASSMASK');
msg = sprintf('Calculation successful!\n\tFull time spent: %3.2f minutes;\n\tTime spent on mask calculation: %3.2f minutes;\n\tMask data size: %3.2f MBytes',...
    (Tm+Tb+Ts+To)/60,(Tm)/60,s.bytes/1024^2);

mailbot('alexander.kharlan@yaliny.com',msg,'MATLAB MAIL BOT', 'RESULT.csv');

clear T Tb Tm Ts To s msg succ ans vtau dtau stau HR a e om Om in u dOm du ddu br_m br_f Vminday DMBmo SGBmo





