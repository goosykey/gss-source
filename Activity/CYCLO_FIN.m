function [ POWERS, SUNANGLES, SUBSAT, Rxyz ] = CYCLO_FIN( EPHEMERIS, jd, Pout, Nreal, map, mapdata )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

R0 = 6378.14; deg = pi/180;
latlim = mapdata{1}; lonlim = mapdata{2};
latstep = mapdata{3}; lonstep = mapdata{4};
M = mapdata{5}; N = mapdata{6}; cliffname = mapdata{7};

lats1 = 90:-1:-90;
lons1 = -180:179;
[latmap1, lonmap1] = ndgrid(lats1,lons1); 

[Msmall, Nsmall] = size(latmap1);

lats = latlim(1):-latstep:latlim(2);
lons = lonlim(1):lonstep:lonlim(2);

map = [zeros(round((90-latlim(1))/latstep),numel(lons)); map; zeros(round((latlim(2)+90)/latstep),numel(lons))];

lats = 90:-latstep:-90;
lons = -180:lonstep:180;
[latmap, lonmap] = ndgrid(lats,lons);

areamap = (R0*pi/180*latstep*lonstep)^2./cosd(latmap);
densmap = map./areamap; densmap(map(:)==0)=0;
densmapsmall = imresize(densmap>750,[Msmall Nsmall],'bilinear');
densmap750small = double(densmapsmall > 4.8e04);
densmap750small(densmap750small==0)=nan;

M = round(N/2);

t = EPHEMERIS(:,1); NN = numel(t); tc = t(2)-t(1);
jt = jd + t/86400;

a = EPHEMERIS(:,2);
h = a-R0;

e = EPHEMERIS(:,3);
om = EPHEMERIS(:,4);
Om = EPHEMERIS(:,5);
in = EPHEMERIS(:,6);
u = EPHEMERIS(:,7);

[Rxyz,~] = randv(a,e,in,Om,om,u-om);
SUBSAT = zeros(NN,2);

%% RECEPTION CONDITIONS

% mapradii = zeros(M,N,3);
% mapradii(:,:,1) = R0*cosd(latmap).*cosd(lonmap);
% mapradii(:,:,2) = R0*cosd(latmap).*sind(lonmap);
% mapradii(:,:,3) = R0*sind(latmap);

%[ ~, ~, ~, ~, Omask, ~, ~, ~, ~] = decompose( BIGASSMASK );

Omaskfull = oceanALL(latmap1,lonmap1,5);
cellno = latmap1;
cellno(:) = 1:numel(cellno);



cliffmap = csvread(cliffname);
cliffmapsmall = imresize(double(cliffmap),[Msmall Nsmall],'nearest');
%cliffmapsmall(cliffmapsmall==0)=nan;
cliffmap = boolean(imresize(double(cliffmap),[M N],'bilinear'));



weat1 = latmap1*0;
[weat2,~] = getweather(weat1,Omaskfull);

global temp;

rcvgood = zeros(NN,1);

fprintf(' -> Processing scanning conditions...\n');
figure;
for i = 1:NN
    if mod(i,round(600/tc))==1 || i==1
        [weat2,weatb] = getweather(weat2,Omaskfull);
        weatbBIG = boolean(imresize(double(weatb),[M N],'nearest'));
        CONDIT = weatbBIG | cliffmap | (densmap>4.8e04) ; % BAD is 1
    end
    
    
    Rxyzecf = R3(gstime(jt(i)))*Rxyz(:,i);
    SUBSAT(i,1) = asin(Rxyzecf(3)/(R0+h(i)))/deg;
    SUBSAT(i,2) = atan2(Rxyzecf(2),Rxyzecf(1))/deg;
    latsat = SUBSAT(i,1); lonsat = SUBSAT(i,2);

    latmapcr = latmap1(abs(latmap1-latsat) < 8.5);
    lonmapcr = lonmap1(abs(latmap1-latsat) < 8.5);
    cellnocr = cellno(abs(latmap1-latsat) < 8.5);

    di = R0*deg*distance(latmapcr,lonmapcr,SUBSAT(i,1),SUBSAT(i,2));
    zoneh = di>750;
    temp = zoneh;
    
    zones = boolean(ones(Msmall,Nsmall));
    zones(cellnocr) = zoneh;
    zones = boolean(imresize(double(zones),[M N],'nearest'));
    CONDITnow = CONDIT | zones; % 1 is BAD
    Sg = sum(areamap(~CONDITnow(:)));
    Sz = sum(areamap(~zones(:)));
    
    rcvgood(i) = Sg/Sz;
    if rcvgood(i) == 1
        rep = 'A';
    else
        rep = sprintf('%g',floor(rcvgood(i)*10));
    end
    
    fprintf(rep);
    
    if mod(i,round(600/tc))==1 || i==1
        weatb = double(weatb);
        weatb(~weatb)=nan;
    end
    %temp = cliffmapsmall;
    surf(lons1,lats1,weatb,'FaceColor','c'); hold on;
    surf(lons1,lats1,cliffmapsmall);
    %surf(lons1,lats1,densmap750small,'FaceColor','m');
    
    plot(lonsat,latsat,'*r','MarkerSize',8);
    
        ttt = 0:1:360;
        gamma = atan(750/R0);
        zg = R0*cos(gamma);
        xg = R0*sin(gamma)*cosd(ttt);
        yg = R0*sin(gamma)*sind(ttt);
        GG = zeros(3,length(ttt));
        GG(1,:) = xg; GG(2,:) = yg; GG(3,:) = zg;
        XYZ = R3(pi/2-SUBSAT(i,2)*pi/180)*R1(pi/2-latsat*pi/180)*GG;
        FIS = asin(XYZ(3,:)/R0)*180/pi;
        LAS = atan2(XYZ(2,:),XYZ(1,:))*180/pi;
        LAS(abs(LAS(:)-lonsat)>180) = LAS(abs(LAS(:)-lonsat)>180) + 360;
        LAS((LAS(:)-lonsat)>180) = LAS((LAS(:)-lonsat)>180) - 720;
        plot(LAS,FIS,'r', 'LineWidth',2);
    
    toptopstring = rep;
    text(0.5,1,toptopstring,'Units','normalized', 'VerticalAlignment','bottom',...
        'HorizontalAlignment','center', 'FontSize',10, 'Color','black');
    view(0,90);
    drawnow;hold off;
    %pause(0.3);

    
end
fprintf('\n');

%temp = SUBSAT;

%% POWER CALCULATION

fprintf(' -> Calculating dissipated power......\n');

SUBSATshift = SUBSAT(2:end,2); SUBSATshift(end+1) = SUBSAT(1,2);

SUBSATshift = SUBSATshift - SUBSAT(:,2);


SUBSATplot = SUBSAT;
SUBSATplot(SUBSATshift>10,:) = nan;

figure
plot(SUBSATplot(:,2),SUBSATplot(:,1),'r','LineWidth',1);
hold on; grid on;
plot_google_map;

Nsubs = sum(Nreal,2);


Pout(Pout>700) = 700;

Pscan = ones(NN,1)*1525;
condit = (Nsubs<3000) .* rcvgood ;
Pscan = Pscan - 1210*condit;


Pout_MAX = [5.8 21 52 92 144 208 283 370 468 578 700];
Pout_EIRP = [0.04 0.76 3.87 12.2 29.8 62 114 195 313 478 700];
Nsa = [1 4 9 16 25 36 49 64 81 100 121];

Pout_EIRP_MATRIX = repmat(Pout_EIRP,NN,1);
Pout_MATRIX = repmat(Pout,1,11);

Pout_COMPARE = Pout_EIRP_MATRIX - Pout_MATRIX;

Pout_COMPARE(Pout_COMPARE<0) = inf;

[~,indexlist] = min(Pout_COMPARE,[],2);

Nsa_list = Nsa(indexlist)';
Pout_MAX_list = Pout_MAX(indexlist)';

Pout1 = Pout;
Pout1(Pout1==0)=1e-7;

Pdis = (Pout1.*(sqrt(Pout_MAX_list./Pout1)/0.26 - 1) + 4.52*Nsa_list + Pscan)*1/0.85;
Pcon = Pdis+Pout;

POWERS = [Pout, Pdis, Pcon];

%% SOLAR CONDITION ANALYSIS

fprintf(' -> Analyzing solar conditions...\n');

au1 = 149597870.700; %km

Rsun = sun(jt)*au1;

SUNANGLES = zeros(NN,3);

for i = 1:NN
    R1sun = Rsun(:,i) - Rxyz(:,i);
    R1suntsw = xyz2tsw(R1sun, Om(i), in(i), u(i));
    SUNANGLES(i,1) = acos(R1suntsw(3) / norm(R1suntsw));
    SUNANGLES(i,2) = atan2(R1suntsw(2),R1suntsw(1));
    lit = sight(Rsun(:,i),Rxyz(:,i),'e');
    if strcmp(lit, 'yes')
        SUNANGLES(i,3)=1;
    else
        SUNANGLES(i,3)=0;
    end
end



end

