%% CONSTANTS

R0 = 6378.137;
deg = pi/180;
au1 = 149598000; %km
mu = 398600.44;


%% INIT
alt = 200; %km


t = 0:10:86400*90;

a = R0+alt; e = 0; om = 0; in = astro.sunsyn(alt)*deg; % specify initial orbit parameters
Om = 0*deg; u = 0*deg;

ini0 = [a e om Om in u];

gd0 = [2018 3 20 12 0 0];

jd0 = astro.jday(gd0(1),gd0(2),gd0(3),gd0(4),gd0(5),gd0(6));
jt = jd0+t/86400;

Nt = numel(t);
dt = t(2) - t(1);

%% IMP

eph = astro.J2pert (ini0,t);

[R,~] = math.randv(eph(:,1),eph(:,2),eph(:,5),eph(:,4),eph(:,3),eph(:,6)-eph(:,3));
ephcart = R';

rsun = astro.sun(jt);

%% SUNLIGHT

rsunkm = rsun' * au1;

lit = 0 * t';

for i = 1:Nt
    litst = astro.sight(ephcart(i,:),rsunkm(i,:),'e');
    if strcmp(litst,'yes')
        lit(i) = 1;
    else
        lit(i) = 0;
    end
end

plot(t/60, lit);
hold on; grid;

aa = diff(find(diff(lit)))*dt;
period = 2* pi* sqrt(a.^3./mu);
aaa = abs(aa(2)-aa(1));
shadow = (period - aaa)/2;
perc = shadow/period;

clear a e i om in u Om ini0 R rsun rsunkm litst dt alt aaa