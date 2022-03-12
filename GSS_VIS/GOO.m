function GOO 

gd0 = [2017 06 04 17 53 45]; % SET UP STARTTIME (UTC)
%gd0 = [2017 03 20 10 28 0]; % SET UP STARTTIME (UTC)

R0 = 6378.137;
deg = pi/180;
elev = 40;
POI = [33,-84];
lat = POI(1);
lon = POI(2);

jd = astro.jday(gd0(1),gd0(2),gd0(3),gd0(4),gd0(5),gd0(6));

cdata = imread('EM_low.jpg');

%% INITIAL DEPLOYMENT
alt = 1150;
elevcrit = elev; % deg
alphacrit = asind(sind(90 + elevcrit) * R0 / (R0 + alt));
beta = 90-elevcrit-alphacrit;

a = R0+alt; e = 0; om = 0; in = 53*deg; % specify initial orbit parameters
Om = 0*deg; u = 0*deg;
ini0 = [a e om Om in u];

% Rcov = R0 * sind(beta);

[ ~, CARTN1, ~ ] = sat.seedconstel( ini0, jd, 32, 50, 360 );



vis_drawstate( CARTN1, jd, 32, 50 );
%add3d.all_coverage( CARTN, alphacrit, ':w', 'c' );

%% FINAL DEPLOYMENT BIG
alt = 1110;
elevcrit = elev; % deg
alphacrit = asind(sind(90 + elevcrit) * R0 / (R0 + alt));
beta = 90-elevcrit-alphacrit;

a = R0+alt; in = 53.8*deg; % specify initial orbit parameters
Om = 5*deg; u = 3*deg;
ini0 = [a e om Om in u];

% Rcov = R0 * sind(beta);

[ ~, CARTN2, ~ ] = sat.seedconstel( ini0, jd, 32, 50, 360 );

vis_drawstate( CARTN2, jd, 32, 50, cdata, 'og',[0.75 0.75 0.75], 'g', false );
%add3d.all_coverage( CARTN, alphacrit, ':w', 'c' );

%% FINAL DEPLOYMENT 1
alt = 1130;
elevcrit = elev; % deg
alphacrit = asind(sind(90 + elevcrit) * R0 / (R0 + alt));
beta = 90-elevcrit-alphacrit;

a = R0+alt; in = 74*deg; % specify initial orbit parameters
Om = 5*deg; u = 3*deg;
ini0 = [a e om Om in u];

% Rcov = R0 * sind(beta);

[ ~, CARTN3, ~ ] = sat.seedconstel( ini0, jd, 8, 50, 360 );

vis_drawstate( CARTN3, jd, 8, 50, cdata, 'oy','y', 'y', false );
%add3d.all_coverage( CARTN, alphacrit, ':w', 'c' );

%% FINAL DEPLOYMENT 2
alt = 1130;
elevcrit = elev; % deg
alphacrit = asind(sind(90 + elevcrit) * R0 / (R0 + alt));
beta = 90-elevcrit-alphacrit;

a = R0+alt; in = 81*deg; % specify initial orbit parameters
Om = 25*deg; u = 4*deg;
ini0 = [a e om Om in u];

% Rcov = R0 * sind(beta);

[ ~, CARTN4, ~ ] = sat.seedconstel( ini0, jd, 5, 75, 360 );

vis_drawstate( CARTN4, jd, 5, 75, cdata, 'om','m', 'm', false );
%add3d.all_coverage( CARTN, alphacrit, ':w', 'c' );

%% FINAL DEPLOYMENT 3
alt = 1325;
elevcrit = elev; % deg
alphacrit = asind(sind(90 + elevcrit) * R0 / (R0 + alt));
beta = 90-elevcrit-alphacrit;

a = R0+alt; in = 70*deg; % specify initial orbit parameters
Om = 45*deg; u = 4*deg;
ini0 = [a e om Om in u];

% Rcov = R0 * sind(beta);

[ ~, CARTN5, ~ ] = sat.seedconstel( ini0, jd, 6, 75, 360 );

vis_drawstate( CARTN5, jd, 6, 75, cdata, 'ob','b', 'b', false );
%add3d.all_coverage( CARTN, alphacrit, ':w', 'c' );


%% GLOB

global bull;

bull = [CARTN1; CARTN2; CARTN3; CARTN4; CARTN5];
bull(:,4:6) = [];
[Nsat, ~] = size(bull);

X = R0 * cosd(lat) * cosd(lon);
Y = R0 * cosd(lat) * sind(lon);
Z = R0 * sind(lat);
XYZ = repmat([X Y Z],[Nsat,1]);
R1s = bull-XYZ;
DOT = dot(XYZ,R1s,2);
XYZabs = sqrt(sum(XYZ.^2,2));
R1sabs = sqrt(sum(R1s.^2,2));
angles = acosd(DOT./XYZabs./R1sabs);
anglecrit = 90-elev;
Ncool = sum(angles < anglecrit)
add3d.all_coverage( bull(angles < anglecrit,:), alphacrit, ':w', 'c' );

%% FINALLY

add3d.majorcities(jd, true);
add3d.sunlight(jd);
add3d.moonpoint(jd);

% dim = [0 0 .3 .3];
% 
% toplstring = sprintf('Diameter: %3.2f km', 2*Rcov);
% annotation('textbox',dim,'String',toplstring,'FitBoxToText','on',...
%     'color', 'w', 'edgecolor', 'w', 'fontsize',7);



end