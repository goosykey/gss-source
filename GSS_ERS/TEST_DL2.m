clear
clc

R0 = 6378.14;
deg = pi/180;

lon = -180:10:180;
lat = -90:5:90;

ttt = zeros(numel(lat),numel(lon));
dlt = zeros(numel(lat),numel(lon));

GS = [69.374366, 88.133709
67.657603, 53.046734
51.406888, 128.142310
77.913685, 16.231174
-51.654133, -69.378787
-21.957492, -43.267948
-33.942142, 18.541134
43.420484, -5.752250
55.717807, 37.442979
54.771077, 20.438429
-31.991908, 116.086949
46.766040, -71.149253
49.381489, -122.892055
49.831284, -97.358232
-5.231855, -80.629570
-34.581975, -56.13440
49.874073, 8.663150
55.971237, 92.850055
62.716823, 40.208079
17.299641, 96.536660
14.746748, -17.466233
35.643995, 51.511861];

Ngs = numel(GS(:,1));

%parpool(2);
%fprintf('Parallel latitude scanning begins... \n');

[~, ~, ~, ~,~, ~, ephcart, gscart] = ...
            access1('time', 0:30:20e04, 'gs', GS); % initialize ephemeris

for q = 1:numel(lon)
    fprintf('Lon = %3.2f deg -> ',lon(q));
    tic
    for i = 1:numel(lat)
        fprintf('|'); %[EVENTTABLE, EVENTGS, MAXNOACCESS, MAXDL, MAXALL, packtable, ephcart, gscart]
        [~, ~, ttt(i,q), ~,dlt(i,q), ~,~,~] = ...
            access1('poi',[lat(i),lon(q)], 'time', 0:30:20e04, 'error', 0,...
            'gs',GS,'gscart', gscart, 'ephcart', ephcart, 'fullaccess',1);
    end
    to = toc;
    fprintf(' - Done in %3.2f minutes.\n',to/60);
end

ttt = ttt/60;
dlt = dlt/60;

visstyle = 2;

switch visstyle
    case 1
        figure;
        contourf(lon,lat,dlt,10);
        colorbar;
        hold on;
        grid;
        load coastlines;
        geoshow(coastlat, coastlon, 'Color', 'white', 'LineWidth', 1.5);
        plot(GS(:,2),GS(:,1),'*r');
    case 2
        figure('Color','white')
        ax = worldmap('world');
        setm(ax,'MLabelParallel',-90)
        setm(ax,'MLabelLocation',90)
        contourm(lat,lon,dlt,10,'Fill','on','LineColor','black');
        load coastlines;
        geoshow(coastlat, coastlon, 'Color', 'white', 'LineWidth', 1.5);
        contourcbar('peer',ax,'Location','southoutside');
        geoshow(GS(:,1),GS(:,2),'color','red','marker','*','displaytype','point');
    case 3
        % cov GS

        elev = 5; alt = 750;
        alphacrit = asind(sind(90 + elev) * R0 / (R0 + alt));
        beta = 90-elev-alphacrit;
        Scov = 2*pi*R0^2*(1-cosd(beta));
        Rcov = sqrt(Scov/pi)*1e3;
        
        cdata = imread('EM_hd.jpg');
        gd0 = [2018 3 20 12 0 0];
        jd = astro.jday(gd0(1),gd0(2),gd0(3),gd0(4),gd0(5),gd0(6));
        [ H ] = vis_earthdraw(jd, 'axes', 'none', 'prime', 'none','quality','hd', 'cdata', cdata);
        Xgs = 1e3*R0 .* cosd(GS(:,1)) .* cos(GS(:,2)*deg - astro.gstime(jd));
        Ygs = 1e3*R0 .* cosd(GS(:,1)) .* sin(GS(:,2)*deg - astro.gstime(jd));
        Zgs = 1e3*R0 .* sind(GS(:,1));
        plot3(Xgs, Ygs, Zgs, 'or','markerfacecolor', 'c');
        for i = 1:Ngs
            cc = [Xgs(i),Ygs(i),Zgs(i)];
            add3d.elem.plotCircle3D(cc,cc,Rcov, '-g')
        end
       
end


% cov GS

% elev = 5; alt = 750;
% alphacrit = asind(sind(90 + elev) * R0 / (R0 + alt));
% beta = 90-elev-alphacrit;
% Scov = 2*pi*R0^2*(1-cosd(beta));
% Rcov = sqrt(Scov/pi)*1e3;
% 
% cdata = imread('EM_hd.jpg');
% gd0 = [2018 3 20 12 0 0];
% jd = astro.jday(gd0(1),gd0(2),gd0(3),gd0(4),gd0(5),gd0(6));
% [ H ] = vis_earthdraw(jd, 'axes', 'none', 'prime', 'none','quality','hd', 'cdata', cdata);
% Xgs = 1e3*R0 .* cosd(GS(:,1)) .* cos(GS(:,2)*deg - astro.gstime(jd));
% Ygs = 1e3*R0 .* cosd(GS(:,1)) .* sin(GS(:,2)*deg - astro.gstime(jd));
% Zgs = 1e3*R0 .* sind(GS(:,1));
% plot3(Xgs, Ygs, Zgs, 'or','markerfacecolor', 'c');
% for i = 1:Ngs
%     cc = [Xgs(i),Ygs(i),Zgs(i)];
%     add3d.elem.plotCircle3D(cc,cc,Rcov, '-g')
% end

%delete(gcp('nocreate'));
