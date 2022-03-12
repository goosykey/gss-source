function GOO_lower 

gd0 = [2017 06 04 17 53 45]; % SET UP STARTTIME (UTC)
%gd0 = [2017 03 20 10 28 0]; % SET UP STARTTIME (UTC)

R0 = 6378.137;
deg = pi/180;
elev = 35;

jd = astro.jday(gd0(1),gd0(2),gd0(3),gd0(4),gd0(5),gd0(6));

cdata = imread('EM_low.jpg');

%% LOWER 1
alt = 345.6; in = 53*deg; num = 2547;

COORDZ = (2*rand(num,1)-1)*sin(in);
HH = sqrt(ones(num,1) - COORDZ.^2);

COORDXY = 2*rand(num,2)-1;
COORDXY = COORDXY ./ repmat(sqrt(sum(COORDXY.^2,2)),[1,2]) .* repmat (HH,[1,2]);

COORD1 = [COORDXY,COORDZ] .* (R0 + alt);

[ H ] = vis_earthdraw(jd, 'axes', 'none', 'prime', 'none','quality','hd', 'cdata', cdata);

X = COORD1(:,1)*1e03;
Y = COORD1(:,2)*1e03;
Z = COORD1(:,3)*1e03;

plot3(X,Y,Z,'or', 'markerfacecolor', [0.75, 0.75, 0.75]); % bull


%% LOWER 2
alt = 340.8; in = 48*deg; num = 2478;

COORDZ = (2*rand(num,1)-1)*sin(in);
HH = sqrt(ones(num,1) - COORDZ.^2);

COORDXY = 2*rand(num,2)-1;
COORDXY = COORDXY ./ repmat(sqrt(sum(COORDXY.^2,2)),[1,2]) .* repmat (HH,[1,2]);

COORD2 = [COORDXY,COORDZ] .* (R0 + alt);

X = COORD2(:,1)*1e03;
Y = COORD2(:,2)*1e03;
Z = COORD2(:,3)*1e03;

plot3(X,Y,Z,'og', 'markerfacecolor', [0.75, 0.75, 0.75]); % bull


%% LOWER 3
alt = 335.9; in = 42*deg; num = 2493;

COORDZ = (2*rand(num,1)-1)*sin(in);
HH = sqrt(ones(num,1) - COORDZ.^2);

COORDXY = 2*rand(num,2)-1;
COORDXY = COORDXY ./ repmat(sqrt(sum(COORDXY.^2,2)),[1,2]) .* repmat (HH,[1,2]);

COORD3 = [COORDXY,COORDZ] .* (R0 + alt);

X = COORD3(:,1)*1e03;
Y = COORD3(:,2)*1e03;
Z = COORD3(:,3)*1e03;

plot3(X,Y,Z,'ob', 'markerfacecolor', [0.75, 0.75, 0.75]); % bull


%% FINALLY

add3d.majorcities(jd, true);
add3d.sunlight(jd);
add3d.moonpoint(jd);

% dim = [0 0 .3 .3];
% 
% toplstring = sprintf('Diameter: %3.2f km', 2*Rcov);
% annotation('textbox',dim,'String',toplstring,'FitBoxToText','on',...
%     'color', 'w', 'edgecolor', 'w', 'fontsize',7);

%% GLOB

global bull;

bull = [COORD1; COORD2; COORD3];


end