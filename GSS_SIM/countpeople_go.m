[ map, mapdata ] = readmap( 'glp15ag15.asc' );    % read population density map
jd0 = astro.jday(2018,1,1,12,0,0);
STEP = 60;

alt = 650;
deg = pi/180;
R0 = 6378.14;

a = R0+alt; e = 0; om = 0; in = astro.sunsyn(alt)*deg; % specify initial orbit parameters
Om = -45*deg; u = 0*deg;

ini0 = [a e om Om in u];
eph0 = astro.J2pert(ini0,43200,STEP);

t = eph0(:,1); lt = length(t);
a = eph0(:,2);
e = eph0(:,3);
om = eph0(:,4);
Om = eph0(:,5);
in = eph0(:,6);
u = eph0(:,7);

jt = jd0 + t/86400;
thetag = astro.gstime(jt);

elevcrit = 10;

alphacrit = asind(sind(90 + elevcrit) * R0 / (R0 + alt));
beta = 90-elevcrit-alphacrit;

Rcov = R0 * sind(beta);
dcrit = R0*sind(beta)/sind(alphacrit); % critical swath distance

latlim = mapdata{1}; lonlim = mapdata{2};
latstep = mapdata{3}; lonstep = mapdata{4};

if latlim(1)>latlim(2)
    latstep = -latstep;
end

LATS = latlim(1):latstep:latlim(2);
LONS = lonlim(1):lonstep:lonlim(2);

[latmap, lonmap] = ndgrid(LATS,LONS);

[M,N] = size(map);

R1ecfmap = R0.*cosd(latmap).*cosd(lonmap);
R2ecfmap = R0.*cosd(latmap).*sind(lonmap);
R3ecfmap = R0.*sind(latmap);

[R, ~] = math.randv(a,e,in,Om,om,u-om);
R = R';
Recf = R * 0;
tete = cell(0);

for i = 1:numel(t)
    if mod(t(i),1800) == 0
        tete = [tete, sprintf('%3.1f',t(i)/3600)];
    end
end

PEOPLE = t*0;


for i = 1:numel(t)
    fprintf('|');
    Recfhere = (math.R3(thetag(i)) * R(i,:)')';
    Recf(i,:) = Recfhere;

    R1ecfsat = repmat(Recfhere(1),M,N);
    R2ecfsat = repmat(Recfhere(2),M,N);
    R3ecfsat = repmat(Recfhere(3),M,N);
    
    DISTMAP = sqrt( (R1ecfsat-R1ecfmap).^2 + (R2ecfsat-R2ecfmap).^2 + (R3ecfsat-R3ecfmap).^2 );
    PEOPLE(i) = sum(sum(map(DISTMAP < dcrit)));
    
    if t(i)/3600 == 1.0
        lath = asind(Recfhere(3)/sqrt(sum(Recfhere.^2)));
        lonh = atan2d(Recfhere(2),Recfhere(1));
        fprintf('\n%g : %g : %g\n',PEOPLE(i), lath, lonh);
        test_dm = DISTMAP;
        boo = test_dm < dcrit;
        maa = map.*boo;
        figure
        surf (LONS, LATS, int8(boo))
    end
end
fprintf('\n');






lateph = asind(Recf(:,3)./sqrt(sum(Recf.^2,2)));
loneph = atan2d(Recf(:,2),Recf(:,1));
deltalon = loneph(2:end) - loneph(1:end-1);
loneph(deltalon > 200) = nan;
lateph(deltalon > 200) = nan;

figure;
plot_google_map;
plot(loneph,lateph);
plot(loneph(mod(t,1800)==0),lateph(mod(t,1800)==0),'or',...
    'markersize',6,'markerfacecolor','red');
text(loneph(mod(t,1800)==0),lateph(mod(t,1800)==0),tete,...
    'verticalalignment','bottom','color','red','fontsize',8);

figure;
plot(t/3600, PEOPLE); 
hold on; grid on;


