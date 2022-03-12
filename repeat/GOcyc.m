%% Moscow orbit - option 1
%  Big camera, weak ADCS

%% INITIALIZE

iniMSC;

%% CONSTANTS

au1 = 149597870.700; %km
sunconst = 1367; %W/m^2
R0 = 6378.137;
deg = pi/180;
mu = 398600.44;


%% IMPORT .XLSX

[~, ~, XL] = xlsread('BOOB.xlsx',3,'A1:AK150');

%% PIECES OF HARDWARE, IMPORT CONSUMPTION DATA

cou = 3;

while (1)
    if isnan( XL{cou,2} )
        break
    end
    cou = cou + 1;
end

Nhard = cou-3; % total number of unique hardware with 2 weirds in the end
HARDcat = cell2mat(XL(3:2+Nhard,1)); % Nhard x 1 double - categories of hardware
HARDno = cell2mat(XL(3:2+Nhard,3)); % Nhard x 1 double - number of pieces
HARDsim = cell2mat(XL(3:2+Nhard,4)); % Nhard x 1 double - number of pieces working simultaneously

CONSon  = cell2mat(XL(3:Nhard,5)); % Nhard-2 x 1 double - consump. ON
CONSon((Nhard-1):(Nhard)) = nan; % 2 nans for shit
CONSsb  = cell2mat(XL(3:2+Nhard,6)); % Nhard x 1 double - consump. StandBy - 2 nans in the end
CONSoff = cell2mat(XL(3:2+Nhard,7)); % Nhard x 1 double - consump. OFF - 2 last nums = effcy
HARDvoltage = cell2mat(XL(3:2+Nhard,8)); % Nhard x 1 double - voltage
HARDvoltage_trans_multipl = cell2mat(XL(3:2+Nhard,9)); % Nhard x 1 double - nnon-bus voltage multiplier

CONSon_total  = HARDsim .* CONSon  .* HARDvoltage_trans_multipl;
CONSon_total(isnan(CONSon_total)) = 0;
CONSsb_total  = HARDsim .* CONSsb  .* HARDvoltage_trans_multipl;
CONSsb_total(isnan(CONSsb_total)) = 0;
CONSoff_total = HARDsim .* CONSoff .* HARDvoltage_trans_multipl;
CONSoff_total(isnan(CONSoff_total)) = 0;

USAGE_coef = cell2mat(XL(3:2+Nhard,14)); % Nhard x 1 double - usage coefficient for partial

CONSpartial_total = (1-USAGE_coef) .* CONSsb_total + USAGE_coef .* CONSon_total;

eff1 = CONSoff(Nhard-1); % power supply components efficiency
eff2 = CONSoff(Nhard);

%% CONSUMPTION MODES

cou = 15;

while (1)
    if isnumeric( XL{1,cou} ) && isnumeric( XL{1,cou+1} ) && isnumeric( XL{1,cou+2} )
        break
    end
    cou = cou + 3;
end

Nmodes = (cou - 15) / 3;

if Nmodes ~= round(Nmodes)
    error ('Consumption modes are fucked up.');
end

% Define mode names:

modeADL = 0;
modeACT = 0;
modePDL = 0;
modeSOL = 0;
modeUMB = 0;

for j = 1:Nmodes
    currenttitle = XL(1,12 + 3*j);
    if contains(currenttitle, 'ADL')
        modeADL = j;
    elseif contains(currenttitle, 'ACT')
        modeACT = j;
    elseif contains(currenttitle, 'PDL')
        modePDL = j;
    elseif contains(currenttitle, 'SOL')
        modeSOL = j;
    elseif contains(currenttitle, 'UMB')
        modeUMB = j;
    end
end

if modeADL == 0 || modeACT == 0 || modePDL == 0 || modeSOL ==0 || modeUMB == 0
    error('Consumption mode titles are shitty.');
end

% Define consumptions in modes:

modeonoff = zeros(Nhard,Nmodes);
modecons  = zeros(Nhard,Nmodes);


for i = 1:Nhard-2
    for j = 1:Nmodes
        cur = XL{i+2,12 + 3*j};
        switch cur
            case 'раб.' % 1
                modeonoff(i,j) = 1;
                modecons (i,j) = CONSon_total(i);
            case 'х.х.' % 2
                modeonoff(i,j) = 2;
                modecons (i,j) = CONSsb_total(i);
            case 'откл.' % 3
                modeonoff(i,j) = 3;
                modecons (i,j) = CONSoff_total(i);
            case 'частичн.' % 4
                modeonoff(i,j) = 4;
                modecons (i,j) = CONSpartial_total(i);
            otherwise
                error ('Some hardware is in fucked up mode.');
        end
    end
    
end

for j = 1:Nmodes % consumption for power supply things
    modeonoff(Nhard-1,j) = 1;
    modeonoff(Nhard,j) = 1;
    modecons(Nhard-1,j) = sum(modecons(1:Nhard-2,j))*(1/eff1-1);
    modecons(Nhard,j)   = sum(modecons(1:Nhard-1,j))*(1/eff2-1);
end



%% POWER SUPPLY SYSTEM PARAMETERS

cou = 2 + Nhard + 1;

while 1
    if ~ isnan(XL{cou,2})
        break
    end
    cou = cou + 1;
end

lifetime  = XL{cou,3}; % years
solarpanels = XL{cou+5,3};
solarpanel_width = XL{cou+6,3};
solarpanel_length = XL{cou+7,3};

solarpanel_area = solarpanel_width * solarpanel_length;
effective_area_cells = XL{cou+9,3};
eff_cells = 0.28*0.8 * (1-0.003)^lifetime;

effective_area = solarpanels * solarpanel_area *...
    effective_area_cells * eff_cells;

adcs_prec_solar = XL{cou+12,3}; % DEG!

Psolar_full = sunconst * effective_area;

eff_charge = XL{cou+25,3};

battery_single_Ah = XL{cou+26,3}; % A*h
battery_single_voltage = XL{cou+27,3}; % Volts

battery_single_cap = battery_single_Ah * battery_single_voltage;

batteries_per_line = XL{cou+29,3}; 

battery_line_voltage = batteries_per_line * battery_single_voltage; % Volts
battery_line_cap = battery_line_voltage * battery_single_Ah; % W*h

batteries_lines = XL{cou+32,3};

battery_full_cap = battery_line_cap * batteries_lines; % W*h


%% INITIAL PARAMETERS

tick = 5; % STEP LENGTH (sec)

transmitter_br = 3e08; % bit/s
camera_br = 2.169e09; % bit/s

radius_circle = 100; % km
elev_crit = 30; % deg

POI = [
    26.8687900000000 100.220720000000 % Lijiang
    45.4642700000000 9.18951000000000 % Milano
    37.2911100000000 127.008890000000 % Suwon-Si
    37.5660000000000 126.978400000000 % Seoul
    37.6563900000000 126.835000000000 % Goyang-Si
    37.4564600000000 126.705150000000 % Incheon
    55.7522200000000 37.6155600000000 % Moscow
    55.1540200000000 61.4291500000000 % Chelyabinsk
    41.8500300000000 -87.6500500000000 % Chicago
    ]; 

POD = [
    %55.75222000 37.61556000 % Moscow
    %59.85000000 29.05000000 % ~St. Petersburg
    %43.114688, 131.956137 % Vladivostok
    %53.349777, 83.769289  % Barnaul
    78.532679, 15.783943 % Spitzbergen
    ];


%% Derived parameters

a = R0 + 560.993907;
Tkepl = 2*pi*sqrt(a^3/mu);

battery_cycles = lifetime*365*86500/Tkepl;
battery_disch_lim = (12345.3/(battery_cycles-395.219)^0.586299 - 1.41611)/100;

nPOI = numel(POI)/2;
nPOD = numel(POD)/2;
POIecf = R0 * [cosd(POI(:,1)).*cosd(POI(:,2)) cosd(POI(:,1)).*sind(POI(:,2)) sind(POI(:,1))];
PODecf = R0 * [cosd(POD(:,1)).*cosd(POD(:,2)) cosd(POD(:,1)).*sind(POD(:,2)) sind(POD(:,1))];

%% GENERATE EPHEMERIS

[ t,y,yfil ] = PEG( jd0, tick, t_fin, 0.5, Y, Y, 'P');

p = y(:,1);
pfil = yfil(:,1);
l1 = y(:,2);
l1fil = yfil(:,2);
l2 = y(:,3);
l2fil = yfil(:,3);
e = sqrt(l1.^2 + l2.^2);
efil = sqrt(l1fil.^2 + l2fil.^2);
a = p./(1-e.^2);
afil = pfil./(1-efil.^2);
om = atan2(l2,l1);
omfil = atan2(l2fil,l1fil);
Om = y(:,4);
Omfil = yfil(:,4);
in = y(:,5);
infil = yfil(:,5);
u = y(:,6);

jt = jd0 + t/86400;

%% R and V, ECI -> ECF

fprintf(' ==> Calculating positions & stuff... \n');

[Rxyz, V] = math.randv(a,e,in,Om,om,u-om);

A = math.R3v(astro.gstime(jt));
n = numel(jt);
Recf = zeros(3,n);


for i = 1:n
    Recf(:,i) = A(:,:,i) * Rxyz(:,i);
end


%% SOLAR CONDITIONS

fprintf(' ==> Calculating solar conditions... \n');

Rsun = astro.sun(jt)*au1;

angletoS = zeros(n,1);
iflit = zeros(n,1);

for i = 1:n
    R1sun = Rsun(:,i) - Rxyz(:,i);
    R1suntsw = math.xyz2tsw(R1sun, Om(i), in(i), u(i));
    angletoS(i) = acos(R1suntsw(2) / norm(R1suntsw));
    lit = astro.sight(Rsun(:,i),Rxyz(:,i),'e');
    if strcmp(lit, 'yes')
        iflit(i) = 1;
    else
        iflit(i) = 0;
    end
end


%% POWER BALANCE

fprintf(' ==> Calculating power balance... \n');

battery_now = battery_full_cap;
dataonboard_now = 0; % bit

CONSUMPTION = zeros(n,1);
PRODUCTION  = zeros(n,1);
battery_log  = zeros(n,1);
dataonboard_log  = zeros(n,1);
modelog = zeros(n,1);

for i = 1:n
    
    % distances to Points of Interest (real and subsat)
    Recfrpt = repmat(Recf(:,i)',[nPOI 1]);
    R1POI = Recfrpt - POIecf;
    DISTtoPOI = sqrt(sum(R1POI.^2, 2));
    Rhere = sqrt(sum(Recf(:,i).^2));
    subsat = [asind(Recf(3,i)/Rhere) atan2d(Recf(2,i),Recf(1,i))];
    SUBSATtoPOI = R0*deg*distance(subsat(1), subsat(2), POI(:,1), POI(:,2));
    
    % distances and angles to Points of Downlink
    Recfrpt = repmat(Recf(:,i)',[nPOD 1]);
    R1POD = Recfrpt - PODecf;
    DISTtoPOD = sqrt(sum(R1POD.^2, 2));
    ANGLEStoPOD = 90 - acosd(dot(PODecf,R1POD,2)/R0./DISTtoPOD); % elev at POD
    
    % DETERMINE THE MODE
    if any(SUBSATtoPOI < radius_circle)
        takingpics = true;
    else
        takingpics = false;
    end
    
    if any(ANGLEStoPOD > elev_crit)
        transmitting_available = true;
    else
        transmitting_available = false;
    end
    
    if dataonboard_now > 0
        transmitting_needed = true;
    else
        transmitting_needed = false;
    end
    
    %main if - choose mode:
    
    if takingpics
        if transmitting_available
            modenow = modeADL;
        else
            modenow = modeACT;
        end
    else % endif if transmitting_available - pics
        if transmitting_available
            if transmitting_needed
                modenow = modePDL;
            else
                if iflit(i)
                    modenow = modeSOL;
                else
                    modenow = modeUMB;
                end
            end
        else
            if iflit(i)
                modenow = modeSOL;
            else
                modenow = modeUMB;
            end
        end % endif if transmitting_available - no pics
    end % endif takingpics
    % mode choose fin
    
    CONSUMPTION(i) = sum(modecons(:,modenow));
    
    switch modenow
        case modeSOL
            PRODUCTION(i) = Psolar_full;
        case modeUMB
            PRODUCTION(i) = 0;
        otherwise
            PRODUCTION(i) = Psolar_full * cos(angletoS(i)) * iflit(i);
            if PRODUCTION(i) < 0
                PRODUCTION(i) = 0;
            end   
    end
    
    
    BALANCE = PRODUCTION(i) - CONSUMPTION(i)/eff_charge; % Positive power income
    
    battery_now = battery_now + BALANCE * tick/3600;
    if(battery_now > battery_full_cap)
        battery_now = battery_full_cap;
    end
    
    dataonboard_now = dataonboard_now + takingpics * camera_br ...
        - transmitting_available * transmitting_needed * transmitter_br;
    if dataonboard_now < 0
        dataonboard_now = 0;
        disp('shet');
    end
    
    battery_log(i) = battery_now;
    dataonboard_log(i) = dataonboard_now;
    modelog(i) = modenow;
end

%% PLOT RESULTS

fprintf(' ==> Enjoy! \n');

epheplot2d( jd0,t,y );

figure;
area (t/3600,PRODUCTION,'facecolor','y'); hold on; grid;
plot (t/3600,CONSUMPTION,'r','linewidth',2)

figure;
area (t/3600,battery_log,'facecolor','r'); hold on; grid;

figure;
area (t/3600,dataonboard_log/8/1024^3,'facecolor','b'); hold on; grid;

% figure
% plot(t,modelog)
% 











