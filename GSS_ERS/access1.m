function [EVENTTABLE, EVENTGS, MAXNOACCESS, MAXDL, MAXALL, packtable, ephcart, gscart] = access1(varargin)
%Returns ground station coordinates in 3D array
%   
%   Initial release: 24.05.2018 by A. Kharlan
%       Returns a 3D array, thus the name xD
%   
%   INPUT:
%       ini0: 1 x 6 - [a e om Om in u] for 1st satellite
%       elev: critical elevation, deg
%       POI : [lat, lon] of POI, deg
%       
%   Output:
%       gscart : 


%% CONSTANTS

R0 = 6378.137;
deg = pi/180;

%% PRE-DEFINE DEFAULT VALUES
alt = 550;
elev = 17;
elev_gs = 5;
pla = 7;
spp = 6;
fangle = 161; % in deg

gd0 = [2018 3 20 12 0 0];

visualize = false;

phase_error = false;

t = 0:30:1e05;

a = R0+alt; e = 0; om = 0; in = astro.sunsyn(alt)*deg; % specify initial orbit parameters
Om = -22.5*deg; u = 0*deg;

ini0 = [a e om Om in u];

POI = [30 0];

GS = [59.875255, 30.313437 % St. Petersburg
    64.735630, 177.495562 % Anadyr
    78.147539, 16.774086 % Svalbard
    ]; 

fullaccess = 0; % no full access
ephmode = 0; % no sat ephemeris defined
gsmode = 0; % no gs ephemeris defined

silentmode = 0;


%% VARARGIN CHECK

nvarargs = length(varargin);

if rem(nvarargs,2) ~= 0
    error('Invalid number of arguments' )
end

while ~isempty(varargin)
    switch lower(varargin{1})
        case 'poi'
            POI = varargin{2};
        case 'gs'
            GS = varargin{2};
        case 'ini0'
            ini0 = varargin{2};
        case 'elev'
            elev = varargin{2};
        case 'elev_gs'
            elev_gs = varargin{2};
        case 'pla'
            pla = varargin{2};
        case 'spp'
            spp = varargin{2};
        case 'time'
            t = varargin{2};
        case 'starttime'
            gd0 = varargin{2};
        case 'visualize'
            visualize = varargin{2};
        case 'error'
            phase_error = varargin{2};
        case 'fullaccess'
            fullaccess = varargin{2};
        case 'ephcart'
            ephcart = varargin{2}; % OVERRIDES PLA & SPP
            ephmode = 1;
        case 'gscart'
            gscart = varargin{2}; % OVERRIDES GS
            gsmode = 1;
        case 'silentmode'
            silentmode = varargin{2};
        otherwise
            error(['Unexpected option: ' varargin{1}])
    end
    varargin(1:2) = [];
end


%% INITIALIZE EVERYTHING ELSE

alphacrit = asind(sind(90 + elev) * R0 / (R0 + alt));

jd0 = astro.jday(gd0(1),gd0(2),gd0(3),gd0(4),gd0(5),gd0(6));
jt = jd0+t/86400;

Nt = numel(t);

if ~ephmode
    Nsat = pla*spp;
else
    Nsat = numel(ephcart(:,1,1));
end

if ~gsmode
    Ngs = numel(GS(:,1));
else
    Ngs = numel(gscart(:,1,1));
end

if phase_error && ephmode
    if ~ silentmode
        warning('Cannot use predefined ephemeris if "phase_error" is set to TRUE. Ephemeris will be ignored.');
    end
    ephmode = 0;
end


%% SEED CONSTELLATION

if ~ephmode % if not using predefined ephemeris
    [ STATEN1, CARTN1, ~ ] = sat.seedconstel( ini0, jd0, pla, spp, fangle );
    
    if visualize
        vis_drawstate( CARTN1, jd0, pla,spp );
        add3d.all_coverage( CARTN1, alphacrit, ':w', 'c' );
    end
    
    if phase_error
        error_vect = (2*rand(numel(STATEN1(:,6)),1) - 1) * 10*deg;
        STATEN1(:,6) = STATEN1(:,6) + error_vect;
    end
end

%% PROPAGATE CONSTELLATION, GS, and POI

if ~ephmode % if not using predefined ephemeris
    [~, ephcart] = astro.J2circ_cube(STATEN1,t); % cartesian ephemeris cube for satellites
end

if ~gsmode % if not using predefined gs
    [gscart] = astro.gs_cube(GS, t, jd0); % cartesian ephemeris cube for GS
end
[POIcart] = astro.gs_cube(POI,t,jd0); % cartesian ephemeris 'cube' for POI, dim1 = 1

%% ACCESS TIMES SATELLITES--POI

ephcart = ephcart(:,:,1:3); % drop velocities
EVENTTABLE = ephcart(:,:,1)*0;

POIcube = repmat(POIcart,[pla*spp, 1, 1]); % comparison cube

r12cube = ephcart - POIcube;

%Rabs = sqrt(sum(ephcart.^2, 3));
DISTANCES = sqrt(sum(r12cube.^2, 3));

GAMMAS = 180-acosd(dot(POIcube, r12cube, 3) ./ (R0.*DISTANCES) );
EVENTTABLE(GAMMAS >= 90+elev) = 1;

%% FIND MAX TIME WITH NO ACCESS

EVENT1 = boolean(sum(EVENTTABLE,1));

ff = find(EVENT1);
maxz = max(diff(ff));

MAXNOACCESS = (t(2)-t(1)) * maxz;

%% ACCESS TIMES SATELLITES--GS[i]

EVENTGS = repmat(EVENTTABLE*0,[1,1,Ngs]);


for i = 1:Ngs
    EVENTGS1 = ephcart(:,:,1)*0;
    
    gs1cube = repmat(gscart(i,:,:),[pla*spp, 1, 1]); % comparison cube
    r12cube = ephcart - gs1cube;
    DISTANCES = sqrt(sum(r12cube.^2, 3));
    
    GAMMAS = 180-acosd(dot(gs1cube, r12cube, 3) ./ (R0.*DISTANCES) );
    EVENTGS1(GAMMAS >= 90+elev_gs) = 1;
    EVENTGS(:,:,i) = EVENTGS1;
end

EVENTGS = sum(EVENTGS,3);
EVENTGS(EVENTGS ~= 0) = 1;

%% FIND MAX DOWNLINK TIME

MAXDL = nan;

if ~fullaccess
    MAXTIMESAT = zeros(1,Nsat);
    
    for i = 1:Nsat % for each satellite
        POIaccess = EVENTTABLE(i,:);
        GSaccess = EVENTGS(i,:);
        LINE = POIaccess*0;
        packagenow = 0;
        for j = 1:Nt % for each point of time
            if POIaccess(j)
                packagenow = 1;
            end
            if GSaccess(j)
                packagenow = 0;
            end
            LINE(j) = packagenow;
        end
        ff = find(~LINE);
        maxz = max(diff(ff));
        MAXTIMESAT(i) = (t(2)-t(1)) * maxz;
    end
    
    MAXDL = max(MAXTIMESAT);
end


%% FULL ACCESS ANALYSIS

MAXALL = nan;
packtable = nan;

if fullaccess
    packtable = zeros(Nt,3); packtable(:,1) = (1:Nt)';
    BIGTABLE = cell(Nsat,Nt);
    poolnow = [];
    for j = 1:Nt % for each point of time / package
        poolnow = [poolnow,j];
        if EVENT1(j) % if sat-POI access here
            satshere = find(EVENTTABLE(:,j)); % find which sats are accessible
            i = satshere(1); % assign package pool to first available satellite
            BIGTABLE{i,j} = [BIGTABLE{i,j} , poolnow]; % assign
            packtable(poolnow,2) = j; % stamp upload time on the packages of the pool
            poolnow = []; % empty the pool
        end
    end
    
    for i = 1:Nsat % for each satellite
        poolnow = [];
        for j = 1:Nt % for each point of time
            poolnow = [poolnow,BIGTABLE{i,j}]; % add packages to satellite's pool
            if EVENTGS(i,j) % if GS is accessible by sat i in time j
                packtable(poolnow,3) = j; % stamp download time on the packages of the pool
                poolnow = []; % empty the pool
            end
        end
    end
    
    maxz = max(packtable(:,3) - packtable(:,1));
    MAXALL = (t(2)-t(1)) * maxz;
    
end








end