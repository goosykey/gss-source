function [TABLE, whichcity] = ad2(varargin)
%AD1 Summary of this function goes here
%   Detailed explanation goes here


global elevs_cube lit_cube elevsun_cube angle_between_cube POI_cube mag_cube jt eph_cube

%% CONSTANTS

R0 = 6378.137;
deg = pi/180;
au1 = 149598000; %km

%% PRE-DEFINE DEFAULT VALUES
alt = 550;
elevcrit = 13;
dt = 15;

gd0 = [2018 2 20 12 0 0];

t = 0:dt:6000;

a = R0+alt; e = 0; om = 0; in = 60*deg;
in = astro.sunsyn(alt)*deg; % specify initial orbit parameters
Om = 60*deg; u = 0*deg;

ini0 = [a e om Om in u];

census = 5e06;

silentmode = false;
dspot = 50;

showpic = 0;

%% VARARGIN CHECK

nvarargs = length(varargin);

if rem(nvarargs,2) ~= 0
    error('Invalid number of arguments' )
end


while ~isempty(varargin)
    switch lower(varargin{1})
        case 'census'
            census = varargin{2};
        case 'dspot'
            dspot = varargin{2};
        case 'ini0'
            ini0 = varargin{2};
        case 'alt'
            alt = varargin{2};
        case 'elevcrit'
            elevcrit = varargin{2};
        case 'time'
            t = varargin{2};
        case 'starttime'
            gd0 = varargin{2};
        case 'silentmode'
            silentmode = varargin{2};
        case 'showpic'
            showpic = varargin{2};
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


%% MODIFY POPULATION - POP IN THE SPOT
try
    load addata.mat poplatlon cdbmodified;
catch ER
    rethrow(ER)
end


if ~exist('poplatlon','var')
    [cdbmodified, poplatlon] = popmodify(census, dspot);
end

cdb = cdbmodified;

Nc = numel(poplatlon(:,1));

%% PROPAGATE, FIND SUN

eph = astro.J2pert (ini0,t);

eph_cube = permute(eph,[3,1,2]);
eph_cube = repmat(eph_cube,[Nc,1,1]);

[R,V] = math.randv(eph(:,1),eph(:,2),eph(:,5),eph(:,4),eph(:,3),eph(:,6)-eph(:,3));
ephcart = R';
vcart = V';

ephcart_cube = permute(ephcart,[3,1,2]);
ephcart_cube = repmat(ephcart_cube,[Nc,1,1]);

POIPOI = poplatlon(:,2:3);

POI_cube = astro.gs_cube(POIPOI,t,jd0);

[rsun, RAsun, DEsun] = astro.sun(jt);

% rsunnorm = [cos(DEsun) .* cos(RAsun)
%             cos(DEsun) .* sin(RAsun)
%             sin(DEsun)];
Dsun = sqrt(sum(rsun.^2));
rsunnorm = rsun ./ repmat(Dsun,[3,1]);

%% ACCESS TIMES SATELLITES--POI

r12_cube = ephcart_cube - POI_cube;

DISTANCES_cube = sqrt(sum(r12_cube.^2,3));

GAMMAS_cube = 180-acosd(dot(POI_cube,r12_cube,3) ./ (R0 .* DISTANCES_cube));

elevs_cube = GAMMAS_cube - 90;

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

lit_cube = repmat(lit',[Nc,1]);

%% SUN @ POI

rsunnorm_cube = permute(rsunnorm,[3,2,1]);
rsunnorm_cube = repmat(rsunnorm_cube,[Nc,1,1]);

POInorm_cube = POI_cube ./ R0;
sunangle_cube = acosd(dot(POInorm_cube,rsunnorm_cube,3));
elevsun_cube = 90-sunangle_cube;

DATETIME = datetime(jt,'convertfrom','juliandate');
dn = datenum(DATETIME);

%% ANGLE BETWEEN

angle_between_cube = acosd(dot(r12_cube,rsunnorm_cube,3)./DISTANCES_cube);

%% OBSERVATION ANGLE

% RADII_cube = sqrt(sum(ephcart_cube.^2,3));
% ANGLES_cube = 90 - acosd(dot(ephcart_cube,r12_cube,3) ./ (RADII_cube.*DISTANCES_cube) );

%% MAGNITUDE

angle_gamma = 180 - angle_between_cube;

mag_cube = GetMagnitude('elev',elevs_cube,'gamma',angle_gamma, 'sm', pi*25);

%angle_gamma(114, 2701:2714)'

%% REVERSED IMAGE THINGY

u_cube = mod(eph_cube(:,:,6),360);


%% LOCAL TIME LIMITATION

timeinpoint = datetime(jt,'convertfrom','juliandate','TimeZone','UTC');
zd = timezone(poplatlon(:,3));

[ZD,TIME] = ndgrid(zd, timeinpoint);

TIME = TIME - ZD/24;

HOURS = TIME.Hour;


%% CHECK CONDITIONS

C1 = elevs_cube > 0;         % Satellite visible from POI
C2 = lit_cube;               % Satellite is lit
C3 = elevsun_cube < 5;       % Sun is below certain elev. @ POI
C4 = angle_between_cube > 25;% Sun doesn't interfere
C5 = mag_cube < -8;          % MAGNITUDE
%C6 = ( (u_cube >= 0) & (u_cube <= 90*deg) ) | u_cube > 270*deg; % reversed image is not shown
C7 = (HOURS >= 6 & HOURS < 11) | (HOURS >= 17 & HOURS < 23);    % local time limitations

SI = C1 & C2 & C3 & C4 & C5 & C7;

% C1(114,2706)
% C2(114,2706)
% C3(114,2706)
% C4(114,2706)
% C5(114,2706)
% C6(114,2706)
% C7(114,2706)

%% CHOOSE CITY FOR THE SHOW

popN = repmat(poplatlon(:,1),[1,Nt]);
SI_checksum = sum(SI,1) > 0;

SIpop = popN.*SI;

[maxpop,maxno] = max(SIpop,[],1);

whichcity = maxno .* SI_checksum;
howmany =  maxpop .* SI_checksum;

%% POINTING ISSUES

CHARTSTRINGS = cell(Nt,3);

TABLE = cell(Nt, 7);
currentshow = 0;

for i = 2:Nt
    if (whichcity(i) ~= whichcity(i-1)) && (whichcity(i) ~= 0)
        currentshow = currentshow+1;
        TABLE{currentshow,1} = currentshow;      
        TABLE{currentshow,7} = i;
        CIT = whichcity(i); currentshowtime = dt;
        satcoord = ephcart(i,:)';
        newcoord = POI_cube(whichcity(i),i,:);
        newcoord = permute(newcoord,[3,1,2]);
        if whichcity(i-1) == 0
            oldcoord = [0;0;0];
        else
            oldcoord = POI_cube(whichcity(i-1),i-1,:);
            oldcoord = permute(oldcoord,[3,1,2]);
        end
        
        r21_old = oldcoord - satcoord;
        r21_new = newcoord - satcoord;
        divbysums = sqrt(sum(r21_old.^2)) .* sqrt(sum(r21_new.^2));
        pointing_angle = acosd(dot(r21_old,r21_new) ./ divbysums);
        
        timehere = TIME(CIT,i);

        CHARTSTRINGS{i,1} = [cdb{CIT,2}, ', ', cdb{CIT,1}];
        CHARTSTRINGS{i,2} = datestr(timehere,'HH:MM');
        CHARTSTRINGS{i,3} = strcat(num2str(round(pointing_angle)), '\circ ');
        
        TABLE(currentshow,2) = CHARTSTRINGS(i,1);
        TABLE{currentshow,3} = datestr(timehere,'dd.mm HH:MM');
        TABLE{currentshow,5} = poplatlon(CIT,1);
        TABLE{currentshow,6} = [num2str(round(pointing_angle)), ' deg'];
        
        for j = (i+1):Nt
            if whichcity(j) == CIT
                currentshowtime = currentshowtime + dt;
            else
                TABLE{currentshow,4} = currentshowtime;
                break
            end
        end
        
    end
end

TABLE(cellfun('isempty',TABLE(:,1)),:) = [];
TABLE(cell2mat(TABLE(:,4))<45,:) = [];

for i = 1:numel(TABLE(:,1))
    TABLE{i,1} = i;
end

TABLE(end+1,:) = {0, 'TOTAL', '', sum(cell2mat(TABLE(:,4))), sum(cell2mat(TABLE(:,5))),nan, nan};

%TABLE(:,7) = [];




%% FIND MEAN

%MEAN = mean(SI.*ANGLES);


%% PLOT RESULTS

fprintf('now plot \n');

if ~silentmode
    DATETIME = datetime(jt,'convertfrom','juliandate');
    dn = datenum(DATETIME);
    plot(DATETIME,howmany);
    hold on;
%     plot(DATETIME,SI.*elevs);
%     xtickformat('dd/MM HH:mm:ss')
    % datetick('x','dd-mmm-yy HH:MM:SS','keepticks','keeplimits');
    grid;
    % ax = gca;
    % ax.XTick = dn(1:100:end);
    xtickangle(60);
    
    xtickformat('MMM dd')
    
    %text(DATETIME,howmany*1.05,CHARTSTRINGS(:,1),'rotation',75,'clipping','on');
    %text(DATETIME,howmany*0.9,CHARTSTRINGS(:,2),'HorizontalAlignment','right','clipping','on');
    %text(DATETIME,howmany*0.9,CHARTSTRINGS(:,3),'rotation',-15,'color','r','clipping','on');
   
end

%% VISUALIZE STATIC

if showpic > 0
    tt = 0:0.01:2*pi;
    q = showpic;
    jd = jt(q);
    c = whichcity(q);
    
    vis_earthdraw(jd,'axes','both','prime','both');
    add3d.sunlight(jd);
    add3d.moonpoint(jd);
    add3d.majorcities(jd,true);
    
    plot3(ephcart_cube(1,q,1)*1e3,ephcart_cube(1,q,2)*1e3,ephcart_cube(1,q,3)*1e3,'or','markerfacecolor','c');
    
    % plot orbit
    NORMAL = cross(ephcart(q,1:3), vcart(q,1:3));
    NORMAL = NORMAL/norm(NORMAL); % normalizing the normal vector
    A = NORMAL(1); B = NORMAL(2); C = NORMAL(3);
    Rsp = norm(ephcart(q,1:3));
    X = Rsp/sqrt(A^2+C^2) * (C*cos(tt) - A*B*sin(tt)/norm(NORMAL));
    Y = Rsp*sqrt(A^2+C^2)*sin(tt)/norm(NORMAL);
    Z = -Rsp/sqrt(A^2+C^2) * (A*cos(tt) + B*C*sin(tt)/norm(NORMAL));
    plot3(X*1e3,Y*1e3,Z*1e3,'r', 'linewidth', 1);
    
    if c ~= 0
        plot3(POI_cube(c,q,1)*1e3,POI_cube(c,q,2)*1e3,POI_cube(c,q,3)*1e3,'or','markerfacecolor','y');
        
        X12 = [POI_cube(c,q,1) POI_cube(c,q,1)+r12_cube(c,q,1)]*1e3;
        Y12 = [POI_cube(c,q,2) POI_cube(c,q,2)+r12_cube(c,q,2)]*1e3;
        Z12 = [POI_cube(c,q,3) POI_cube(c,q,3)+r12_cube(c,q,3)]*1e3;
        plot3(X12,Y12,Z12,'y','LineWidth',2);
        
        XS = [ephcart_cube(c,q,1) ephcart_cube(c,q,1)+rsunnorm_cube(c,q,1)*1e3]*1e3;
        YS = [ephcart_cube(c,q,2) ephcart_cube(c,q,2)+rsunnorm_cube(c,q,2)*1e3]*1e3;
        ZS = [ephcart_cube(c,q,3) ephcart_cube(c,q,3)+rsunnorm_cube(c,q,3)*1e3]*1e3;
        plot3(XS,YS,ZS,'y','LineWidth',2);
        
        add3d.elem.city3d(jd, poplatlon(c,2),poplatlon(c,3),...
            strcat(cdb{c,2}, ', ', cdb{c,1}), 'color', 'y', 'markersize',7,'fontsize',10);
        
        text(1e3*ephcart(q,1),1e3*ephcart(q,2),1e3*ephcart(q,3),...
            strcat('Mag',sprintf(' %3.1f', mag_cube(c,q))),...
            'VerticalAlignment','bottom',...
        'HorizontalAlignment','center', 'color','c', 'FontSize',10);
    end

end

end

