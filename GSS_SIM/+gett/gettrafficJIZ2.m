function [ data1k, data2k, data3k, JIZ1, JIZ2, JIZ3, hexa_c_ll_here ] = gettrafficJIZ2( lat, lon, R, dt, map1, map2, map3, mapdata, netdata )
%GET TRAFFIC USAGE IN ZONE
%     INPUT:
%         VARCOORD: zone definition. Options:
%               1 - by city name  : VARCOORD = {City, R, cdb}
%               2 - circular area : VARCOORD = [latcenter, loncenter, Radius]
%               3 - rectangular zone : VARCOORD = [lat1 lon1 lat2 lon2]
%         VARTIME : time definition.
%               IF two numbers: period is abs(VARTIME(1)) to abs(VARTIME(2))
%               IF one number : VARTIME is initial time, hours. Period is 1 second.
%               for UTC time use NEGATIVE values. For LOCAL time - POSITIVE
%         map	: population distribution map (not density! i.e. 'PEOPLE', NOT 'PEOPLE/sq.km')
%         mapdata : map data obtained thru 'readmap.m';
%         BIGASSMASK : usage intensity masks cell array. See 'intmasknew.m'
%         pieflag : plot zonal pie diagrams (binary flag addition)
%               1 - MAP
%               2 -  all population
%               4 - all sub
%               8 - 2 types of sub
%              16 - voice + dataN + dataS usage
%              32 - map stripped of neighbouring coverage (32 = 33)
%     OUTPUT:
%         vdata, ddata, sdata : data usage
%         SESTABLE : session coordinates and time


if numel (lat) > 1
    latneigh = lat(2:end);
    lat = lat(1);
else
    latneigh = [];
end

if numel (lon) > 1
    lonneigh = lon(2:end);
    lon = lon(1);
else
    lonneigh = [];
end

map = map1 + map2 + map3;

if dt <0
    skipgenerate = true;
    dt = -dt;
else
    skipgenerate = false;
end

%% GET DAILY DISTRIBUTIONS

tau1 = netdata{4};
tau2 = netdata{5};
tau3 = netdata{6};
averate1 = netdata{7};
averate2 = netdata{8};
averate3 = netdata{9};
bitrate = netdata{10};

%% INITIALIZE MAP PARAMETERS

latlim = mapdata{1}; lonlim = mapdata{2};
latstep = mapdata{3}; lonstep = mapdata{4};

lats = latlim(1):-latstep:latlim(2);
lons = lonlim(1):lonstep:lonlim(2);

[latmap, lonmap] = ndgrid(lats,lons);

%% INTEGRATE TIME MASK

tint1 = averate1/86400 * dt;
tint2 = averate2/86400 * dt;
tint3 = averate3/86400 * dt;

%% RANGES NET (polar) or RECT NET (rectangular), COUNT SUBSCRIBERS



maphere = map;

RANGES = 6378.14*pi/180*distance(latmap,lonmap,lat,lon);
maphere(RANGES(:)>=R) = 0;
map1(RANGES(:)>=R) = 0; map2(RANGES(:)>=R) = 0; map3(RANGES(:)>=R) = 0;

subs1 = sum(sum(map1));
subs2 = sum(sum(map2));
subs3 = sum(sum(map3));

subs0 = sum(sum(maphere));

map(abs(lats-lat) > 10, :) = [ ]; map(:, abs(lons-lon) > 10) = [ ];
map1(abs(lats-lat) > 10, :) = [ ]; map1(:, abs(lons-lon) > 10) = [ ];
map2(abs(lats-lat) > 10, :) = [ ]; map2(:, abs(lons-lon) > 10) = [ ];
map3(abs(lats-lat) > 10, :) = [ ]; map3(:, abs(lons-lon) > 10) = [ ];
maphere(abs(lats-lat) > 10, :) = [ ]; maphere(:, abs(lons-lon) > 10) = [ ];
latmap(abs(lats-lat) > 10, :) = [ ]; latmap(:, abs(lons-lon) > 10) = [ ];
lonmap(abs(lats-lat) > 10, :) = [ ]; lonmap(:, abs(lons-lon) > 10) = [ ];

lons (abs(lons-lon) > 10) = [];
lats (abs(lats-lat) > 10) = [];

%% ∆»«Ќя ќ¬ў»Ќј

[ hexacx, hexacy, hexasx, hexasy ] = sx.hexnew( lon, lat, R, 14 );

R0 = 6378.137;
alt = 800;


% найти кеплеровы координаты  ј - соседей - латлон соседей - обрезать
% €чейки - см gettrafficsat

dists1 = distance(lat, lon, hexacy, hexacx)*pi/180*6378.14;

if ~isempty(latneigh)
    for q = 1:numel(latneigh)
        dists = distance(latneigh(q), lonneigh(q), hexacy, hexacx)*pi/180*6378.14;
        hexacx(dists < dists1) = [];
        hexacy(dists < dists1) = [];
        hexasx(dists < dists1, :) = [];
        hexasy(dists < dists1, :) = [];
        dists1(dists < dists1) = [];
    end
    
end

hexa_c_ll_here = [ hexacx, hexacy];

Nhex = numel(hexacx);

subs1k = zeros(Nhex,2);
subs2k = zeros(Nhex,2);
subs3k = zeros(Nhex,2);
subs0k = zeros(Nhex,2);
cou = 0;

fprintf('Go cycle JIZ, total = %g\n',Nhex);
for k = 1:Nhex
    fprintf('|');
    cou = cou+1;
    
    xv = hexasx (k,:);
    yv = hexasy (k,:);
    maphexmask = inpolygon (lonmap,latmap,xv,yv);
    
    subs1k(cou,1) = sum(sum(map1.*maphexmask)); subs1k(cou,2) = k;
    subs2k(cou,1) = sum(sum(map2.*maphexmask)); subs2k(cou,2) = k;
    subs3k(cou,1) = sum(sum(map3.*maphexmask)); subs3k(cou,2) = k;
    
    subs0k(cou,1) = sum(sum(maphere.*maphexmask)); subs0k(cou,2) = k;
end
fprintf('\n');


%% COUNT AND SHOW USAGE DATA

data1 = tint1*subs1;
data2 = tint2*subs2;
data3 = tint3*subs3;

if not(skipgenerate)
    
    fprintf('In area: TOTAL (type1 + type2 + type3) \n\tSubs : %g thsnd (%g + %g + %g)\n',...
        subs0/1000, subs1/1000, subs2/1000, subs3/1000);
    
    fprintf('TOTAL USAGE \n\tType1 : %3.2f KB \n\tType2 : %3.2f KB \n\tType3 : %3.2f KB \n',...
        data1/1024, data2/1024, data3/1024);
    
end

%% START SESSION MODELLING - find MEAN ACTIVE SUBS

if not(skipgenerate)
    
    bitratebytes = bitrate/8;
    
    W1 = data1/bitratebytes;
    W2 = data2/bitratebytes;
    W3 = data3/bitratebytes;
    
    Nactmean1 = W1/dt;
    Nactmean2 = W2/dt;
    Nactmean3 = W3/dt;
    
end

%% ∆изн€ков 2

data1k = tint1.*subs1k(:,1);
data2k = tint2.*subs2k(:,1);
data3k = tint3.*subs3k(:,1);
W1k = data1k./bitratebytes;
W2k = data2k./bitratebytes;
W3k = data3k./bitratebytes;
Nactmean1k = W1k./dt;
Nactmean2k = W2k./dt;
Nactmean3k = W3k./dt;



%fprintf('Average number of active subscribers \n\tVoice : %3.2f\n\tDataN : %3.2f\n\tDataS : %3.2f\n', vNactmean, dNactmean, sNactmean);

%% FIND NUMBER OF SESSIONS BEGINNING

if not(skipgenerate)
    
    Qzone1 = Nactmean1/tau1*dt;
    Qzone2 = Nactmean2/tau2*dt;
    Qzone3 = Nactmean3/tau3*dt;
    
    fprintf('Avg number of active subs & Sesesions begun during specified period \n\tType1 : \t%8.2f\t\t||\t\t%3.2f\n\tType2 : \t%8.2f\t\t||\t\t%3.2f\n\tType3 : \t%8.2f\t\t||\t\t%3.2f\n',...
        Nactmean1, Qzone1, Nactmean2,Qzone2, Nactmean3,Qzone3);
    %fprintf('Sessions begun during specified period \n\tVoice : %3.2f\n\tDataN : %3.2f\n\tDataS : %3.2f\n', vQzone, dQzone, sQzone);
    
end

%% ∆изн€ков 3

Qzone1k = Nactmean1k/tau1*dt;
Qzone2k = Nactmean2k/tau2*dt;
Qzone3k = Nactmean3k/tau3*dt;

JIZ1 = [Qzone1k, subs1k];
JIZ2 = [Qzone2k, subs2k];
JIZ3 = [Qzone3k, subs3k];


if not(skipgenerate)
    
    %% GENERATE THE SESSIONS
    
    % smart roundation
    Qzone1 = floor(Qzone1) + (mod(Qzone1,1)>rand);
    Qzone2 = floor(Qzone2) + (mod(Qzone2,1)>rand);
    Qzone3 = floor(Qzone3) + (mod(Qzone3,1)>rand);
    
    
    SEST1 = zeros(Qzone1,2);
    
    SEST2 = zeros(Qzone2,2);
    
    SEST3 = zeros(Qzone3,2);
    
    
    
    
    
    figure('Name','Sessions beginning during period','NumberTitle','on'); hold on;
    
    
    fprintf('Generating Type1 Sessions... \n');
    for q = 1:Qzone1
        [lonsub,latsub] = math.pinky(lons,lats,map1,3);
        SEST1(q,:) = [latsub,lonsub];
    end
    
    fprintf('Generating Type2 Sessions... \n');
    for q = 1:Qzone2
        [lonsub,latsub] = math.pinky(lons,lats,map2,3);
        SEST2(q,:) = [latsub,lonsub];
    end
    
    fprintf('Generating Type3 Sessions... \n');
    for q = 1:Qzone3
        [lonsub,latsub] = math.pinky(lons,lats,map3,3);
        SEST3(q,:) = [latsub,lonsub];
    end
    
    
    plot(SEST1(:,2),SEST1(:,1),'o', 'MarkerSize', 5, 'MarkerFaceColor', [1 0.5 0.5]);
    plot(SEST2(:,2),SEST2(:,1),'o', 'MarkerSize', 5, 'MarkerFaceColor', [0.25 0.75 1]);
    plot(SEST3(:,2),SEST3(:,1),'o', 'MarkerSize', 5, 'MarkerFaceColor', [0.5 1 0.5]);
    grid;
    
    
    
    R0 = 6378.14;
    t = 0:1:360;
    gamma = atan(R/R0);
    zg = R0*cos(gamma);
    xg = R0*sin(gamma)*cosd(t);
    yg = R0*sin(gamma)*sind(t);
    GG = zeros(3,length(t));
    GG(1,:) = xg; GG(2,:) = yg; GG(3,:) = zg;
    XYZ = math.R3(pi/2-lon*pi/180)*math.R1(pi/2-lat*pi/180)*GG;
    FIS = asin(XYZ(3,:)/R0)*180/pi;
    LAS = atan2(XYZ(2,:),XYZ(1,:))*180/pi;
    LAS(abs(LAS(:)-lon)>180) = LAS(abs(LAS(:)-lon)>180) + 360;
    LAS((LAS(:)-lon)>180) = LAS((LAS(:)-lon)>180) - 720;
    plot(LAS,FIS,'r', 'LineWidth',2);
    
    
    bottomlstring = sprintf('Tick length : %02.2f sec.',dt);
    toprstring = sprintf('Sessions begun now:\n%2.0f T1\n%2.0f T1\n%2.0f T3', Qzone1, Qzone2, Qzone3);
    toplstring = sprintf('Active subs now:\n%2.0f\n%2.0f\n%2.0f', Nactmean1, Nactmean2, Nactmean3);
    
    toptopstring = sprintf('Area subs count: %3.0f. T1: %3.0f T2: %3.0f T3: %3.0f', subs0, subs2, subs2, subs3);
    
    text(1,1,toprstring,'Units','normalized', 'VerticalAlignment','top',...
        'HorizontalAlignment','right', 'FontSize',10, 'Color','red');
    text(0,0,bottomlstring,'Units','normalized', 'VerticalAlignment','bottom',...
        'HorizontalAlignment','left', 'FontSize',10, 'Color','red');
    text(0,1,toplstring,'Units','normalized', 'VerticalAlignment','top',...
        'HorizontalAlignment','left', 'FontSize',10, 'Color','red');
    
    text(0.5,1,toptopstring,'Units','normalized', 'VerticalAlignment','bottom',...
        'HorizontalAlignment','center', 'FontSize',10, 'Color','black');
    
    plot_google_map;
    
    sx.plotem(hexasx, hexasy);
    
end

end

