function TEST1 

gd0 = [2018 03 20 18 0 0]; % SET UP STARTTIME (UTC)
%gd0 = [2017 03 20 10 28 0]; % SET UP STARTTIME (UTC)

R0 = 6378.137;
deg = pi/180;
elev = 15;
POI = [33,-84];
lat = POI(1);
lon = POI(2);

jd = astro.jday(gd0(1),gd0(2),gd0(3),gd0(4),gd0(5),gd0(6));

cdata = imread('EM_low.jpg');

%% INITIAL DEPLOYMENT
alt = 750;
elevcrit = elev; % deg
alphacrit = asind(sind(90 + elevcrit) * R0 / (R0 + alt));
beta = 90-elevcrit-alphacrit;

a = R0+alt; e = 0; om = 0; in = astro.sunsyn(alt)*deg; % specify initial orbit parameters
Om = -22.5*deg; u = 0*deg;
ini0 = [a e om Om in u];

% Rcov = R0 * sind(beta);

[ ~, CARTN1, ~ ] = sat.seedconstel( ini0, jd, 6, 8, 170 );



vis_drawstate( CARTN1, jd, 6, 8 );
add3d.all_coverage( CARTN1, alphacrit, ':w', 'c' );


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