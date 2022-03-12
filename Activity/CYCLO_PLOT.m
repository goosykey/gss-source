function CYCLO_PLOT( gd0, EPHEMERIS, Ntot, Ndes,  Nass, Nreal, POWERS, SUNANGLES )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


jd = jday(gd0(1),gd0(2),gd0(3),gd0(4),gd0(5),gd0(6)); % find starttime Julian Date

R0 = 6378.14; deg = pi/180;
t = EPHEMERIS(:,1); NN = numel(t);
jt = jd + t/86400;
a = EPHEMERIS(:,2);
h = a-R0;
e = EPHEMERIS(:,3);
om = EPHEMERIS(:,4);
Om = EPHEMERIS(:,5);
in = EPHEMERIS(:,6);
u = EPHEMERIS(:,7);
[Rxyz,~] = randv(a,e,in,Om,om,u-om); % ECI cartesian position (skip velocity)
SUBSAT = zeros(NN,2);
for i = 1:NN
    Rxyzecf = R3(gstime(jt(i)))*Rxyz(:,i);
    SUBSAT(i,1) = asin(Rxyzecf(3)/(R0+h(i)))/deg; % subsat latitude
    SUBSAT(i,2) = atan2(Rxyzecf(2),Rxyzecf(1))/deg; % subsat longitude
end

SUBSATshift = SUBSAT(2:end,2); SUBSATshift(end+1) = SUBSAT(1,2);

SUBSATshift = SUBSATshift - SUBSAT(:,2);

SUBSAT(SUBSATshift>10,:) = nan;

tdatetime = t/3600;

%% PLOT GROUND TRACK

figure
plot(SUBSAT(:,2),SUBSAT(:,1),'r','LineWidth',1);
hold on; grid on;
plot_google_map;

%% PLOT SUBS

figure;
area(tdatetime,sum(Ndes,2),'FaceColor','red'); % red area - desired subs
grid on; hold on;
plot(tdatetime,sum(Ntot,2),'-r', 'LineWidth',2); % red dashed - total subs in area
Ndro = Ndes-Nass;
area(tdatetime,sum(Ndes,2) - Ndro(:,2),'FaceColor',[1 0.5 0]); % orange area - desired subs minus dropped dataN subs
area(tdatetime,sum(Nass,2),'FaceColor','yellow'); % yellow area - assisted subs
area(tdatetime,sum(Nreal,2),'FaceColor',[0.8 1 0]); % greenish area - real subs
%datetick('x','dd.mm HH:MM:SS', 'keepticks');

%% PLOT POWERS
figure;
area(tdatetime,SUNANGLES(:,3)*1.2*max(max(POWERS)),'FaceColor','yellow'); % yellow area - IF_SUN
grid on; hold on;
area(tdatetime,POWERS(:,3),'FaceColor',[0 0.5 1]); % blue area - consumed power
area(tdatetime,POWERS(:,2),'FaceColor','cyan'); % cyan area - dissipated power
area(tdatetime,POWERS(:,1),'FaceColor','magenta'); % magenta area - irradiated power
%datetick('x','dd.mm HH:MM:SS', 'keepticks');

end

