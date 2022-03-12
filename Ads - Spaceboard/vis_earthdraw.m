function [ fig, globe, globemisc] = vis_earthdraw(jd, varargin)
%% Textured 3D Earth example
%
% GOOSE SATELLITE SYSTEMS
% 8 Sep 2004 (Ryan Gray)
% Revised 9 March 2006, 31 Jan 2006, 16 Oct 2013

%% VARARGIN CHECKS

%default
axopt = 'ecf';
eqopt = true;
pmopt = 'ecf';

mapq = 'hd';

%image_file = 'http://upload.wikimedia.org/wikipedia/commons/thumb/c/cd/Land_ocean_ice_2048.jpg/1024px-Land_ocean_ice_2048.jpg';
image_file = ['EM_' mapq '.jpg'];

cdata = [];

nvarargs = length(varargin);

if rem(nvarargs,2) ~= 0
    error('Invalid number of arguments' )
end

while ~isempty(varargin)
    switch lower(varargin{1})
        case 'axes'
            axopt = varargin{2};
        case 'equator'
            eqopt = varargin{2};
            try
                eqopt = boolean(eqopt);
            catch ME
                %error('"Equator" expected to be boolean or booleanable' );
                rethrow(ME);
            end
        case 'prime'
            pmopt = varargin{2};
        case 'quality'
            mapq = varargin{2};
            image_file = ['EM_' mapq '.jpg'];
        case 'cdata'
            cdata = varargin{2};
        otherwise
            error(['Unexpected option: ' varargin{1}])
    end
    varargin(1:2) = [];
end

if isempty(cdata)
    cdata = imread(image_file); % Load Earth image for texture map
end

%% Options

space_color = 'k';
npanels = 180;   % Number of globe panels around the equator deg/panel = 360/npanels
alpha   = 1; % globe transparency level, 1 = opaque, through 0 = invisible
%GMST0 = []; % Don't set up rotatable globe (ECEF)
GMST0 = astro.gstime (jd); % Set up a rotatable globe at J2000.0

% Earth texture image
% Anything imread() will handle, but needs to be a 2:1 unprojected globe image.

% Mean spherical earth

erad    = 6378.137; % equatorial radius (KM)
prad    = 6371.0087714; % polar radius (KM)
%erot    = 7.2921158553e-5; % earth rotation rate (radians/sec)

%% Create figure

scrsz = get(groot,'ScreenSize');

fig = figure;
set(gcf, 'Renderer', 'OpenGL', 'OuterPosition',[100,100,scrsz(3)-100,scrsz(4)-100])
%set(H,'Renderer','zbuffer');
vw = get(gca,'View');
set(gca, 'Clipping', 'off'); % do not crop 3d objects
setAxes3DPanAndZoomStyle(zoom(gca),gca,'camera'); % old camera behaviour

if vw(2) == 90
    vw(2)=30;    
end
 cla(fig);
set(fig, 'Color', space_color, 'Name', 'GSS 3D visualization');




hold on;

% Turn off the normal axes

set(gca, 'NextPlot','add', 'Visible','off');

axis equal;

% Set initial view

view(vw);
set(gca, 'XLimMode', 'manual', 'YLimMode', 'manual', 'ZLimMode', 'manual');
set(gca, 'XLim', [-9e+03 9e+03], 'YLim', [-9e+03 9e+03], 'ZLim', [-9e+03 9e+03]);

axis vis3d;

axis([-5e3 5e3 -5e3 5e3 -5e3 5e3])
%set(gca,'dataaspectratio',[1 1 1]);

%% Create wireframe globe

% Create a 3D meshgrid of the sphere points using the ellipsoid function

[x, y, z] = ellipsoid(0, 0, 0, erad, erad, prad, npanels);

globe = surf(x, y, -z, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1], ...
    'SpecularStrength',0.4);

if ~isempty(GMST0)
    hgx = hgtransform;
    set(hgx,'Matrix', makehgtform('zrotate',GMST0));
    set(globe,'Parent',hgx);
end

%% Texturemap the globe

% Set image as color data (cdata) property, and set face color to indicate
% a texturemap, which Matlab expects to be in cdata. Turn off the mesh edges.

set(globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', alpha, 'EdgeColor', 'none');

%% plot misc

% axe frames
if strcmp(axopt, 'both')
    plotaxes('eci',GMST0);
    ecfmisc = plotaxes('ecf',GMST0);
else
    ecfmisc = plotaxes(axopt,GMST0);
end

%equator
if eqopt
    plotequator;
end

%prime
if strcmp(pmopt, 'both')
    plotprime('eci',GMST0);
    prime = plotprime('ecf',GMST0);
else
    prime = plotprime(pmopt,GMST0);
end

globemisc = [ecfmisc;prime];



end



function [misc] = plotaxes(which,th)
%PLOT AXES

switch which
    case 'eci'
        Mx = [0 8e3 ; 0 0 ; 0 0 ]';
        My = [0 0 ; 0 8e3 ; 0 0 ]';
        Mz = [0 0 ; 0 0 ; 0 8e3 ]';
        Tx = [8e3 0 0]; Ty = [0 8e3 0]; Tz = [0 0 8e3];
        tline = {'X_E_C_I', 'Y_E_C_I', 'Z_E_C_I'};
        eciax = plot3(Mx, My, Mz, 'r', 'linewidth',2);
        ecitext = text(Tx, Ty, Tz, tline, 'color','r');
        misc = [eciax; ecitext];
    case 'ecf'
        Mx = [0 8e3*cos(th) ; 0 -8e3*sin(th) ; 0 0 ]';
        My = [0 8e3*sin(th) ; 0 8e3*cos(th) ; 0 0 ]';
        Mz = [0 0 ; 0 0 ; 0 8e3 ]';
        Tx = [8e3*cos(th) -8e3*sin(th) 0];
        Ty = [8e3*sin(th) 8e3*cos(th) 0]; 
        Tz = [0 0 8e3];
        tline = {'X_E_C_F', 'Y_E_C_F', 'Z_E_C_F'};
        ecfax = plot3(Mx, My, Mz, 'g', 'linewidth',2);
        ecftext = text(Tx, Ty, Tz, tline, 'color','g');
        misc = [ecfax; ecftext];
end

end



function plotequator
%PLOT EQUATOR

R0 = 6378.137;
t = 0:0.002:2*pi;

xeq = R0*cos(t);
yeq = R0*sin(t);
zeq = 0*t;

plot3(xeq, yeq, zeq, 'w', 'linewidth',2);
end

function [prime] = plotprime(which,th)
%PLOT PRIME MERIDIAN

R0 = 6378.137;
t = 0:0.002:2*pi;

switch which
    case 'eci'
        xpr = R0*cos(t);
        ypr = 0*t;
        zpr = R0*sin(t);
        prime = plot3(xpr, ypr, zpr, 'w', 'linewidth',2);
    case 'ecf'
        xgr = R0/sin(th) * (sin(th)*cos(th)*sin(t));
        ygr = R0*sin(th)*sin(t);
        zgr = -R0/sin(th) * (sin(th)*cos(t));
        prime = plot3(xgr, ygr, zgr, 'g', 'linewidth',2);
        
end

end
