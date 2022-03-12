function [ H, latsat, area ] = vis_regions1( CART510, jd, cdata )
%UNTITLED4 Summary of this function goes here
%   ALL CALCULATIONS IN KM, PLOT IN METERS!!!

global temp temp1;

temp1 = [];

%% INITIALIZE

R0 = 6378.137; %km
rzone = 675.02;
deg = pi/180;



[ H ] = vis_drawstate( CART510, jd, 17, 30, cdata );
set(H,'units','normalized','position',[.1 .1 .6 .6]);

% 241 is the shit (-> 211 240 242 270 271 300)

x1 = CART510(383,1);
y1 = CART510(383,2);
z1 = CART510(383,3);

R1s = sqrt(x1^2+y1^2+z1^2);
x1s = x1*R0/R1s; y1s = y1*R0/R1s; z1s = z1*R0/R1s; % SUBSAT 1
lat1 = asin(z1/R1s)/deg;
lon1 = atan2(y1,x1)/deg; % THETA SUBSAT 1

EAST1 = [-sind(lon1), cosd(lon1), 0];
NADIR1 = [x1s y1s z1s] - [x1 y1 z1];

plot3(x1*1e3, y1*1e3, z1*1e3, 'pg', 'markersize', 10);

%[349 348 384 382 419 418] boob

COORD6 = CART510([349 384 419 418 382 348],1:3);

apexbetha = []; apextheta = [];

for i = [349 384 419 418 382 348]
    x2 = CART510(i,1);
    y2 = CART510(i,2);
    z2 = CART510(i,3);
    plot3(x2*1e3, y2*1e3, z2*1e3, 'or', 'markerfacecolor', 'g'); % redraw sat
    
    A = 2*(x2-x1); B = 2*(y2-y1); C = 2*(z2-z1);
    ABSABC = sqrt(A^2+B^2+C^2);
    
    t = 0:0.002:2*pi;
    
    x = R0/sqrt(A^2+C^2) * (C*cos(t) - A*B*sin(t)/ABSABC);
    y = R0*sqrt(A^2+C^2)/ABSABC * sin(t);
    z = -R0/sqrt(A^2+C^2) * (A*cos(t) + B*C*sin(t)/ABSABC);
    
    DISTANCESTO6 = zeros(6,length(x));
    for j = 1:6
        DISTANCESTO6(j,:) = sqrt((x-COORD6(j,1)).^2 + (y-COORD6(j,2)).^2 + (z-COORD6(j,3)).^2);
    end
    
    
    DISTANCES = sqrt((x-x2).^2 + (y-y2).^2 + (z-z2).^2); % dist to the current neigh
    x (min(DISTANCESTO6,[],1) < DISTANCES) = [];
    y (min(DISTANCESTO6,[],1) < DISTANCES) = [];
    z (min(DISTANCESTO6,[],1) < DISTANCES) = [];
    
    bethas = asind(z/R0);
    thetas = atan2d(y,x);
    
    %% plot coverage
    R2s = sqrt(x2^2+y2^2+z2^2);
    t = 0:1:360;
    lat2 = asin(z2/R2s)/deg;
    lon2 = atan2(y2,x2)/deg;
    gamma = atan(rzone/R0); 
    zg = R0*cos(gamma);
    xg = R0*sin(gamma)*cosd(t);
    yg = R0*sin(gamma)*sind(t);
    GG = zeros(3,length(t));
    GG(1,:) = xg; GG(2,:) = yg; GG(3,:) = zg;
    XYZ = math.R3(pi/2-lon2*pi/180)*math.R1(pi/2-lat2*pi/180)*GG;
    plot3(XYZ(1,:)*1e3,XYZ(2,:)*1e3,XYZ(3,:)*1e3,'w','LineStyle','--');
     
    
    %% skip if hexagon -> pentagon -> rectangle
    
    if isempty(bethas)
        continue
    end
    
    %% correct arc begin/end
    
    if distance(bethas(1),thetas(1),bethas(end),thetas(end))*deg*R0 < 10 % shitty arc correction
        fprintf(' > Performing shitty arc correction on sat %g\n',i);
        % look for latitude break
        bethasshift = [bethas(end),bethas];
        bethasshift(end) = [];
        bethasdelta = bethas - bethasshift;
        [~,n] = max(abs(bethasdelta),[],2);
        l = length(bethas);
        neworder = [n:l,1:n-1];
        x = x(neworder);
        y = y(neworder);
        z = z(neworder);
        bethas = bethas(neworder);
        thetas = thetas(neworder);
        
    end
    
    %% plot that shit
    
    apexbetha = [apexbetha;bethas(end)];
    apextheta = [apextheta;thetas(end)];
    
    plot3(x*1e3, y*1e3, z*1e3, 'c', 'linewidth',2);
    
    %endstr = sprintf('end%g',i);
    %begstr = sprintf('beg%g',i);
    %text(x(end)*1e3,y(end)*1e3,z(end)*1e3,endstr, 'color','m');
    %text(x(1)*1e3,y(1)*1e3,z(1)*1e3,begstr, 'color','y');
    
    %% bullshit
    
    r1 = [x(end) y(end) z(end)] - [x1 y1 z1];
    r2 = [x(end) y(end) z(end)] - [x1s y1s z1s];
    r3 = [x1s y1s z1s];
    r4 = [x(end) y(end) z(end)];
    
    %PSI = acos(dot(NADIR1,r1)/norm(NADIR1)/norm(r1))/deg;
    PSI = acos(dot(r3,r4)/norm(r3)/norm(r4))/deg;
    PHI = acos(dot(EAST1,r2)/norm(EAST1)/norm(r2))/deg;
    if z(end) < z1s
        PHI = 360-PHI;
    end
    
    temp1 = [temp1;[PSI,PHI]];
    
end

% apextheta = mod(apextheta,360);

% xxx = 1e3*R0*cosd(apexbetha).*cosd(apextheta);
% yyy = 1e3*R0*cosd(apexbetha).*sind(apextheta);
% zzz = 1e3*R0*sind(apexbetha);
% 
% text(xxx,yyy,zzz,{'1';'2';'3';'4';'5';'6'}, 'color','w');

area = areaint(apexbetha,apextheta,wgs84Ellipsoid)/1e6;%*4*pi*R0^2;


centerstr = sprintf('Area = %3.2f sq.km',area);
text(x1s*1e3,y1s*1e3,z1s*1e3,centerstr, 'color','y');
latsat = asin(z1/sqrt(x1^2+y1^2+z1^2));

view([x1,y1,z1]);

%% plot coverage for central sat

t = 0:1:360;
gamma = atan(rzone/R0);
zg = R0*cos(gamma);
xg = R0*sin(gamma)*cosd(t);
yg = R0*sin(gamma)*sind(t);
GG = zeros(3,length(t));
GG(1,:) = xg; GG(2,:) = yg; GG(3,:) = zg;
XYZ = math.R3(pi/2-lon1*pi/180)*math.R1(pi/2-lat1*pi/180)*GG;
plot3(XYZ(1,:)*1e3,XYZ(2,:)*1e3,XYZ(3,:)*1e3,'r','linewidth',2);

%% bullshit

temp = [[lat1,lon1];[apexbetha, apextheta]];

end

