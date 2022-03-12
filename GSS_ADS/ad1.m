function [SI, elevs, lit, elevsun, MEAN] = ad1(varargin)
%AD1 Summary of this function goes here
%   Detailed explanation goes here

%% CONSTANTS

R0 = 6378.137;
deg = pi/180;
au1 = 149598000; %km

%% PRE-DEFINE DEFAULT VALUES
alt = 400;
elevcrit = 13;
xsicrit = 160;

gd0 = [2018 1 20 12 0 0];

t = 0:15:86400 * 60;

a = R0+alt; e = 0; om = 0; in = 48*deg;
%in = astro.sunsyn(alt)*deg; % specify initial orbit parameters
Om = 0*deg; u = 0*deg;

ini0 = [a e om Om in u];

POI = [45 -93];

silentmode = false;

%% VARARGIN CHECK

nvarargs = length(varargin);

if rem(nvarargs,2) ~= 0
    error('Invalid number of arguments' )
end


while ~isempty(varargin)
    switch lower(varargin{1})
        case 'poi'
            POI = varargin{2};
        case 'ini0'
            ini0 = varargin{2};
        case 'alt'
            alt = varargin{2};
        case 'in'
            in = varargin{2};
        case 'elevcrit'
            elevcrit = varargin{2};
        case 'xsicrit'
            xsicrit = varargin{2};
        case 'time'
            t = varargin{2};
        case 'starttime'
            gd0 = varargin{2};
        case 'silentmode'
            silentmode = varargin{2};
        otherwise
            error(['Unexpected option: ' varargin{1}])
    end
    varargin(1:2) = [];
end

%% INITIALIZE EVERYTHING ELSE

alphacrit = asind(sind(90 + elevcrit) * R0 / (R0 + alt));

jd0 = astro.jday(gd0(1),gd0(2),gd0(3),gd0(4),gd0(5),gd0(6));
jt = jd0+t/86400;

Nt = numel(t);

ini0 = [a e om Om in u];


%% PROPAGATE, FIND SUN

eph = astro.J2pert (ini0,t);

[R,~] = math.randv(eph(:,1),eph(:,2),eph(:,5),eph(:,4),eph(:,3),eph(:,6)-eph(:,3));
ephcart = R';

[POIcart] = astro.gs_cube(POI,t,jd0);
POIcart = permute(POIcart,[2,3,1]);

rsun = astro.sun(jt);
Dsun = sqrt(sum(rsun.^2));
rsunnorm = rsun ./ repmat(Dsun,[3,1]);

%% ACCESS TIMES SATELLITES--POI

r12 = ephcart - POIcart;
DISTANCES = sqrt(sum(r12.^2,2));

GAMMAS = 180-acosd(dot(POIcart,r12,2) ./ (R0 .* DISTANCES));
elevs = GAMMAS - 90;

%% SUNLIGHT

rsunkm = rsun' * au1;

lit = 0 * t';

for i = 1:Nt
    litst = astro.sight(ephcart(i,:),rsunkm(i,:),'e');
    if strcmp(litst,'yes')
        lit(i) = 1;
    else
        lit(i) = 0;
    end
end

%% SUN @ POI

POInorm = POIcart ./ R0;
sunangle = acosd(dot(POInorm,rsunnorm',2));
elevsun = 90-sunangle;

%% ANGLE BETWEEN

angle_between = acosd(dot(r12,rsunnorm',2)./DISTANCES);

%% CHECK CONDITIONS

C1 = elevs > 0;         % Satellite visible from POI
C2 = lit;               % Satellite is lit
C3 = elevsun < 5;       % Sun is below certain elev. @ POI
C4 = angle_between > 25;% Sun doesn't interfere

SI = C1 & C2 & C3 & C4;

%% OBSERVATION ANGLE

RADII = sqrt(sum(ephcart.^2,2));
ANGLES = 90 - acosd(dot(ephcart,r12,2) ./ (RADII.*DISTANCES) );

%% FINDMEAN

MEAN = mean(SI.*ANGLES);


%% PLOT RESULTS

if ~silentmode
    DATETIME = datetime(jt,'convertfrom','juliandate');
    dn = datenum(DATETIME);
    plot(DATETIME,SI.*ANGLES);
    hold on;
    plot(DATETIME,SI.*elevs);
    xtickformat('dd/MM HH:mm:ss')
    % datetick('x','dd-mmm-yy HH:MM:SS','keepticks','keeplimits');
    grid;
    % ax = gca;
    % ax.XTick = dn(1:100:end);
    xtickangle(60);
    
    xtickformat('MMM dd')
end

end

