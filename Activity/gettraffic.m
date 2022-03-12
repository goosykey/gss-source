function [ vdata, ddata, sdata, SESTABLE ] = gettraffic( VARCOORD, VARTIME, map, mapdata, BIGASSMASK, netdata, pieflag )
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

%% ARG MANIPULATIONS


if nargin < 7
    pieflag = 0;
end

if pieflag > 63
    pieflag = 63;
end

if iscell(VARCOORD)
    if iscell(VARCOORD{3}) && isnumeric(VARCOORD{2}) && ischar(VARCOORD{1})
        cdb = VARCOORD{3};
        citylist = cdb(:,2);
        ism = ismember(citylist, VARCOORD(1));
        popul = cell2mat(cdb(:,3));
        [~,maxi] = max(popul.*ism);
        lat = cell2mat(cdb(maxi,4));
        lon = cell2mat(cdb(maxi,5));
        polar = true;
        R = VARCOORD{2};
    else
        error('Argument assignment error: VARCOORD must be either: \n\t {City, R, cdb} \t where cdb is a cell array city database;\n\t [lat, lon, R]\n\t [lat1, lon1, lat2, lon2]\n');
    end
elseif isnumeric(VARCOORD) && (numel(VARCOORD)==3 || numel(VARCOORD)==4)
    if numel(VARCOORD)==3
        polar = true;
        lat = VARCOORD(1);
        lon = VARCOORD(2);
        R = VARCOORD(3);
    else
        polar = false;
        fi1 = VARCOORD(1);
        fi2 = VARCOORD(3);
        la1 = VARCOORD(2);
        la2 = VARCOORD(4);
    end
else
    error('Argument type invalid');
end

if numel(VARTIME) >= 2
    if VARTIME(1) * VARTIME(2) < 0 && numel(VARTIME) > 2
        error('Invalid time input.');
    end
end


if VARTIME(1) < 0 % THEN UTC
    if polar
        z = -timezone(lon);
    else
        z = -timezone((la2+la1)/2);
    end
    dt1 = -VARTIME(1) + z;
    if numel(VARTIME) == 2
        dt2 = -VARTIME(2) + z;
    else
        dt2 = dt1 + 1/3600;
    end
else
    dt1 = VARTIME(1);
    if numel(VARTIME) == 2
        dt2 = VARTIME(2);
    else
        dt2 = dt1 + 1/3600;
    end
end


%% DECOMPOSE BIG MASK

[ Vmask, Dmask, Smask, ~, ~, Pmask, SPmask, ~, Zmask] = decompose( BIGASSMASK );

%% GET DAILY DISTRIBUTIONS

vtau = netdata{4};
dtau = netdata{5};
stau = netdata{6};
Vminday = netdata{7};
DMBmo = netdata{8};
SGBmo = netdata{9};
bitrate_m = netdata{10};
bitrate_f = netdata{11};

dt = (dt2-dt1)*3600;

% NORMAL voice daily distrib: SECONDS - SECONDS
Vdayactraw = [5 2.5 1.8 1.8 3 6 9 7.5 8.5 10 9.5 7 5]; 
Vnorm = Vminday*60/trapz(0:2*3600:24*3600,Vdayactraw);
Vdayact = Vdayactraw*Vnorm;
vqq = interp1(0:2:24,Vdayact,0:1/3600:24,'spline');

% NORMAL traffic daily distrib: BYTES - SECONDS
Ddayactraw = [4 1 1 2 4 6 10 10 9 10 9 8 4]; 
Dnorm = DMBmo/30*1024^2/trapz(0:2*3600:24*3600,Ddayactraw);
Ddayact = Ddayactraw*Dnorm;
dqq = interp1(0:2:24,Ddayact,0:1/3600:24,'spline');

% SPEEDY traffic daily distrib: BYTES - SECONDS
Sdayactraw = [8 6 3 2 2.5 3 4 5 6 8 10 10 8];
Snorm = SGBmo/30*1024^3/trapz(0:2*3600:24*3600,Sdayactraw);
Sdayact = Sdayactraw*Snorm;
sqq = interp1(0:2:24,Sdayact,0:1/3600:24,'spline');

%% INITIALIZE MAP PARAMETERS

latlim = mapdata{1}; lonlim = mapdata{2};
latstep = mapdata{3}; lonstep = mapdata{4};

lats = latlim(1):-latstep:latlim(2);
lons = lonlim(1):lonstep:lonlim(2);

[M,N] = size(map);

[latmap, lonmap] = ndgrid(lats,lons);

%% INTEGRATE TIME MASK

vtint = trapz(vqq(round(dt1*3600+1):round(dt2*3600+1)));
dtint = trapz(dqq(round(dt1*3600+1):round(dt2*3600+1)));
stint = trapz(sqq(round(dt1*3600+1):round(dt2*3600+1)));

%% FIT UI MASKS TO MAP SIZE & GET MODIFIED MAPS

Vmask = imresize(double(Vmask)/100,[M N],'bilinear');
Dmask = imresize(double(Dmask)/100,[M N],'bilinear');
Smask = imresize(double(Smask)/100,[M N],'bilinear');

Vmap = map.*Vmask.*Pmask;
Dmap = map.*Dmask.*Pmask;
Smap = map.*Smask.*SPmask;

%% SHAPE SUBSCRIBER MAP

Nsubmap = map.*Pmask;
Ssubmap = map.*SPmask;

%% RANGES NET (polar) or RECT NET (rectangular), COUNT SUBSCRIBERS

maphere = map;

if polar
    RANGES = 6378.14*pi/180*distance(latmap,lonmap,lat,lon);
    maphere(RANGES(:)>=R) = 0;
    Nsubmap(RANGES(:)>=R) = 0; Ssubmap(RANGES(:)>=R) = 0;
    Vmap(RANGES(:)>=R)=0; Dmap(RANGES(:)>=R)=0; Smap(RANGES(:)>=R)=0;
        
else
    INRECT = belongstorect( latmap, lonmap, fi1, fi2, la1, la2 );
    maphere(~INRECT(:)) = 0;
    Nsubmap(~INRECT(:)) = 0; Ssubmap(~INRECT(:)) = 0;
    Vmap(~INRECT(:))=0; Dmap(~INRECT(:))=0; Smap(~INRECT(:))=0;

end

Nsubscr = sum(sum(Nsubmap));
Ssubscr = sum(sum(Ssubmap));

Vmodsubs = sum(sum(Vmap));
Dmodsubs = sum(sum(Dmap));
Smodsubs = sum(sum(Smap));

VmodsubsZ = [sum(Vmap(Zmask(:)==0)) sum(Vmap(Zmask(:)==1)) sum(Vmap(Zmask(:)==2)) sum(Vmap(Zmask(:) >2))];

people = sum(sum(maphere));
people0 = round(sum(maphere(Zmask(:)==0))); Nsubscr0 = round(sum(Nsubmap(Zmask(:)==0))); Ssubscr0 = round(sum(Ssubmap(Zmask(:)==0)));
people1 = round(sum(maphere(Zmask(:)==1))); Nsubscr1 = round(sum(Nsubmap(Zmask(:)==1))); Ssubscr1 = round(sum(Ssubmap(Zmask(:)==1)));
people2 = round(sum(maphere(Zmask(:)==2))); Nsubscr2 = round(sum(Nsubmap(Zmask(:)==2))); Ssubscr2 = round(sum(Ssubmap(Zmask(:)==2)));
people3 = round(sum(maphere(Zmask(:) >2))); Nsubscr3 = round(sum(Nsubmap(Zmask(:) >2))); Ssubscr3 = round(sum(Ssubmap(Zmask(:) >2)));

%% COUNT AND SHOW USAGE DATA

fprintf('In area: TOTAL (DISTANT + RURAL + URBAN + OCEAN) \n\tPopulation : %g thsnd (%g + %g + %g + %g)\n\tNorm subs  : %g (%g + %g + %g + %g)\n\tSpeed subs : %g (%g + %g + %g + %g)\n',...
    people/1000, people0/1000, people1/1000, people2/1000, people3/1000, Nsubscr, Nsubscr0, Nsubscr1, Nsubscr2, Nsubscr3, Ssubscr, Ssubscr0, Ssubscr1, Ssubscr2, Ssubscr3);
fprintf('Average intensity \n\tVoice : %3.2f\n\tDataN : %3.2f\n\tDataS : %3.2f\n', Vmodsubs/Nsubscr, Dmodsubs/Nsubscr, Smodsubs/Ssubscr);


vdata = vtint*Vmodsubs;
ddata = dtint*Dmodsubs;
sdata = stint*Smodsubs;

vdatavg = vdata/Nsubscr/60; % minutes per subsc
ddatavg = ddata/Nsubscr/1024^2; % Mbytes per subsc
sdatavg = sdata/Ssubscr/1024^2; % Mbytes per subsc
fprintf('TOTAL USAGE \n\tVoice : %3.2f hours (%3.3f minutes per sub)\n\tDataN : %3.2f GBytes (%3.3f MBytes per sub)\n\tDataS : %3.2f GBytes (%3.3f MBytes per sub)\n',...
    vdata/3600,vdatavg, ddata/1024^3,ddatavg, sdata/1024^3,sdatavg);

%% START SESSION MODELLING - find MEAN ACTIVE SUBS

dspeed = bitrate_m / 8; % byte/sec
sspeed = bitrate_f / 8; % byte/sec

vW = vdata;
dW = ddata/dspeed;
sW = sdata/sspeed;

vNactmean = vW/dt;
dNactmean = dW/dt;
sNactmean = sW/dt;

%fprintf('Average number of active subscribers \n\tVoice : %3.2f\n\tDataN : %3.2f\n\tDataS : %3.2f\n', vNactmean, dNactmean, sNactmean);

%% FIND NUMBER OF SESSIONS BEGINNING

vQzone = vNactmean/vtau*dt;
dQzone = dNactmean/dtau*dt;
sQzone = sNactmean/stau*dt;

fprintf('Avg number of active subs & Sesesions begun during specified period \n\tVoice : \t%3.2f\t\t||\t\t%3.2f\n\tDataN : \t%3.2f\t\t||\t\t%3.2f\n\tDataS : \t%3.2f\t\t||\t\t%3.2f\n',...
    vNactmean, vQzone, dNactmean,dQzone, sNactmean,sQzone);
%fprintf('Sessions begun during specified period \n\tVoice : %3.2f\n\tDataN : %3.2f\n\tDataS : %3.2f\n', vQzone, dQzone, sQzone);

%% PLOT PIE DIAGRAMS

pies = de2bi(pieflag,6);

if pies(2)
    figure('Name','POPULATION ZONAL DISTRIBUTION','NumberTitle','on'); hold on;
    vec = [people0, people1, people2, people3];
    labels = {'Distant','Rural','Urban','Ocean'};
    explode = [1 0 0 1];
    pie (vec,explode,labels);
    axis off
    hold off;
end

if pies(3)
    figure('Name','SUBSCRIBER ZONAL DISTRIBUTION','NumberTitle','on'); hold on;
    
    Z = [Nsubscr0+Ssubscr0, Nsubscr1+Ssubscr1, Nsubscr2+Ssubscr2, Nsubscr3+Ssubscr3];
    labels = {'Distant','Rural','Urban','Ocean'};
    ax1 = subplot(1,1,1);
    explode = [0 0 0 1];
    pie(ax1,Z,explode,labels)
    title(ax1,'All subs');
    axis off
    hold off;
end

if pies(4)
    figure('Name','SUBSCRIBER TYPE ZONAL DISTRIBUTION','NumberTitle','on'); hold on;
    
    X = [Nsubscr0, Nsubscr1, Nsubscr2, Nsubscr3];
    labels = {'Distant','Rural','Urban','Ocean'};
    ax1 = subplot(1,2,1);
    pie(ax1,X,labels)
    title(ax1,'Normal subs');
    
    Y = [Ssubscr0, Ssubscr1, Ssubscr2, Ssubscr3];
    labels = {'Distant','Rural','Urban','Ocean'};
    ax1 = subplot(1,2,2);
    pie(ax1,Y,labels)
    title(ax1,'Speed subs');
    axis off
    hold off;
end

if pies(5);
    figure('Name','VOICE USAGE ZONAL DISTRIBUTION','NumberTitle','on'); hold on;
    vec = VmodsubsZ;
    explode = [0 0 1 1];
    labels = {'Distant','Rural','Urban','Ocean'};
    pie (vec,explode,labels);
    axis off
    hold off;
end


%% GENERATE THE SESSIONS

% smart roundation
vQzone = floor(vQzone) + (mod(vQzone,1)>rand);
dQzone = floor(dQzone) + (mod(dQzone,1)>rand);
sQzone = floor(sQzone) + (mod(sQzone,1)>rand);

vSEST1 = zeros(vQzone,4);
vSEST2 = zeros(vQzone,2);
dSEST1 = zeros(dQzone,4);
dSEST2 = zeros(dQzone,2);

if isnan(sQzone)
    sQzone = 0;
end
sSEST1 = zeros(sQzone,4);
sSEST2 = zeros(sQzone,2);

vqqz = vqq(round(dt1*3600+1:dt2*3600+1));
dqqz = dqq(round(dt1*3600+1:dt2*3600+1));
sqqz = sqq(round(dt1*3600+1:dt2*3600+1));

vqq_cum = cumsum(vqqz);
dqq_cum = cumsum(dqqz);
sqq_cum = cumsum(sqqz);
lcum = length(vqq_cum);

if pies(1)
    figure('Name','Sessions beginning during period','NumberTitle','on'); hold on;
end

fprintf('Generating Voice Sessions... \n');
for q = 1:vQzone
    counter = q;
    sestype = 1;
    ran = rand*vqq_cum(lcum);
    [~, minicum] = min(abs(vqq_cum-ran));
    minicum = dt1*3600+minicum;
    %temp = vqq_cum-ran;
    hr = idivide(minicum,int16(3600),'floor');
    minute = idivide(mod(minicum,3600),int16(60),'floor');
    second = mod(minicum,60);
    [lonsub,latsub] = pinky(lons,lats,Vmap,3);
    vSEST1(counter,:) = [sestype, hr, minute, second];
    vSEST2(counter,:) = [latsub,lonsub];
end

fprintf('Generating DataN Sessions... \n');
for q = 1:dQzone
    counter = q;
    sestype = 2;
    ran = rand*dqq_cum(lcum);
    [~, minicum] = min(abs(dqq_cum-ran));
    minicum = dt1*3600+minicum;
    hr = idivide(minicum,int16(3600),'floor');
    minute = idivide(mod(minicum,3600),int16(60),'floor');
    second = mod(minicum,60);
    [lonsub,latsub] = pinky(lons,lats,Dmap,3);
    dSEST1(counter,:) = [sestype, hr, minute, second];
    dSEST2(counter,:) = [latsub,lonsub];
end

fprintf('Generating DataS Sessions... \n');
for q = 1:sQzone
    counter = q;
    sestype = 3;
    ran = rand*sqq_cum(lcum);
    [~, minicum] = min(abs(sqq_cum-ran));
    minicum = dt1*3600+minicum;
    hr = idivide(minicum,int16(3600),'floor');
    minute = idivide(mod(minicum,3600),int16(60),'floor');
    second = mod(minicum,60);
    [lonsub,latsub] = pinky(lons,lats,Smap,3);
    sSEST1(counter,:) = [sestype, hr, minute, second];
    sSEST2(counter,:) = [latsub,lonsub];
end

SESTABLE = [[double(vSEST1),vSEST2];[double(dSEST1),dSEST2];[double(sSEST1),sSEST2]];

if pies(1) || pies(6) % PLOT IF FLAG SPECIFIED
    plot(vSEST2(:,2),vSEST2(:,1),'o', 'MarkerSize', 5, 'MarkerFaceColor', [1 0.5 0.5]);
    plot(dSEST2(:,2),dSEST2(:,1),'o', 'MarkerSize', 5, 'MarkerFaceColor', [0.25 0.75 1]);
    plot(sSEST2(:,2),sSEST2(:,1),'o', 'MarkerSize', 5, 'MarkerFaceColor', [0.5 1 0.5]);
    grid;
    
    
    if polar
        R0 = 6378.14;
        t = 0:1:360;
        gamma = atan(R/R0);
        zg = R0*cos(gamma);
        xg = R0*sin(gamma)*cosd(t);
        yg = R0*sin(gamma)*sind(t);
        GG = zeros(3,length(t));
        GG(1,:) = xg; GG(2,:) = yg; GG(3,:) = zg;
        XYZ = R3(pi/2-lon*pi/180)*R1(pi/2-lat*pi/180)*GG;
        FIS = asin(XYZ(3,:)/R0)*180/pi;
        LAS = atan2(XYZ(2,:),XYZ(1,:))*180/pi;
        LAS(abs(LAS(:)-lon)>180) = LAS(abs(LAS(:)-lon)>180) + 360;
        LAS((LAS(:)-lon)>180) = LAS((LAS(:)-lon)>180) - 720;
        plot(LAS,FIS,'r', 'LineWidth',2);
    else
        plot([lon(1) lon(1) lon(2) lon(2) lon(1)],[lat(1) lat(2) lat(2) lat(1) lat(1)],'r', 'LineWidth',2);
    end
    
    if polar && ~pies(6)
        [ plotcover ] = getcoverage( lat, lon, 600, R, netdata, true );
        for i = 2:7
            FIS = plotcover(:,i,1);
            LAS = plotcover(:,i,2);
            plot(LAS,FIS,'b', 'LineWidth',2, 'LineStyle','--');
        end
    end
    
    timehr = floor(dt1);
    timemin = mod(dt1,1)*60;
    
    bottomrstring = sprintf('%g:%g:%02.2f local time', timehr, floor(timemin), mod(timemin,1)*60);
    bottomlstring = sprintf('Tick length : %02.2f sec.',(dt2-dt1)*3600);
    toprstring = sprintf('Sessions begun now:\n%2.0f V\n%2.0f D\n%2.0f S', vQzone, dQzone, sQzone);
    toplstring = sprintf('Active subs now:\n%2.0f\n%2.0f\n%2.0f', vNactmean, dNactmean, sNactmean);
    
    toptopstring = sprintf('Area population: %3.0f. Subs: %3.0f N / %3.0f S', people, Nsubscr, Ssubscr);
    
    text(1,0,bottomrstring,'Units','normalized', 'VerticalAlignment','bottom',...
        'HorizontalAlignment','right', 'FontSize',10, 'Color','red');
    text(1,1,toprstring,'Units','normalized', 'VerticalAlignment','top',...
        'HorizontalAlignment','right', 'FontSize',10, 'Color','red');
    text(0,0,bottomlstring,'Units','normalized', 'VerticalAlignment','bottom',...
        'HorizontalAlignment','left', 'FontSize',10, 'Color','red');
    text(0,1,toplstring,'Units','normalized', 'VerticalAlignment','top',...
        'HorizontalAlignment','left', 'FontSize',10, 'Color','red');
    
    text(0.5,1,toptopstring,'Units','normalized', 'VerticalAlignment','bottom',...
        'HorizontalAlignment','center', 'FontSize',10, 'Color','black');
    
    plot_google_map;
    
end

end

