%% CONST

R0 = 6378.14;
deg = pi/180;
mu = 398600.44; % km^3/s^2
au1 = 150e06;

alt = 600; % km

%% INI

dt = 30;
t = 0:dt:86400*30; % 1 MONTH MODELLING TIME

jd0 = astro.jday(2017, 3, 20, 10, 28, 0);
jt = jd0 + t/86400;

INI_STATE = [R0+alt, 0, 0, 0, astro.sunsyn(alt)*deg 0*deg];

T = 2*pi*sqrt((R0+alt)^3/mu);

%% ORBIT PROPAGATION

[ eph ] = astro.j2circ( INI_STATE, t );
eph = permute(eph,[2,3,1]);

[r,v] = math.randv(eph(:,1), eph(:,2), eph(:,5), eph(:,4), eph(:,3), eph(:,6)); % Keplerian -> Cartesian
r = r'; v = v';

%% SUN

[rsun, rtascsun, declsun] = astro.sun(jt); 
rsun = rsun';
rsun1 = rsun;
rsun = rsun ./ repmat(sqrt(sum(rsun.^2,2)),[1,3]);

[ F81,F107, Kp ] = astro.getsolar( jt, 3 );

