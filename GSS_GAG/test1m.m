%% CONSTANTS

R0 = 6378.137;
deg = pi/180;
au1 = 149598000; %km
mu = 398600.44;

%% PRE-DEFINE DEFAULT VALUES
lifetime = 5;

alt = 600;
elevcrit = 10;
dt = 15;

gd0 = [2018 3 20 12 0 0];

t = 0:dt:86400 * 30;

a = R0+alt; e = 0; om = 0;
in = astro.sunsyn(alt)*deg; % specify initial orbit parameters
Om = 15*deg; u = 0*deg;

T = 2*pi*sqrt(a^3/mu);
v = sqrt(mu/a);

ini0 = [a e om Om in u];

GS = [55 37];

Nc = numel(GS(:,1));


%% INITIALIZE EVERYTHING ELSE

alphacrit = asind(sind(90 + elevcrit) * R0 / (R0 + alt));

jd0 = astro.jday(gd0(1),gd0(2),gd0(3),gd0(4),gd0(5),gd0(6));
jt = jd0+t/86400;

Nt = numel(t);

%% PROPAGATE, FIND SUN

eph = astro.J2pert (ini0,t);

eph_cube = permute(eph,[3,1,2]);
eph_cube = repmat(eph_cube,[Nc,1,1]);

[R,V] = math.randv(eph(:,1),eph(:,2),eph(:,5),eph(:,4),eph(:,3),eph(:,6)-eph(:,3));
ephcart = R';
vcart = V';

ephcart_cube = permute(ephcart,[3,1,2]);
ephcart_cube = repmat(ephcart_cube,[Nc,1,1]);

GS_cube = astro.gs_cube(GS,t,jd0);

[rsun, RAsun, DEsun] = astro.sun(jt);
Dsun = sqrt(sum(rsun.^2));
rsunnorm = rsun ./ repmat(Dsun,[3,1]);

%% ACCESS TIMES SATELLITES--GS

r12_cube = ephcart_cube - GS_cube;

DISTANCES_cube = sqrt(sum(r12_cube.^2,3));

GAMMAS_cube = 180-acosd(dot(GS_cube,r12_cube,3) ./ (R0 .* DISTANCES_cube));

elevs_cube = GAMMAS_cube - 90;

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

lit_cube = repmat(lit',[Nc,1]);

%% GAGEL

SUNLIGHT = sum(lit_cube(1,1:round(T/dt)),2) * dt;
SHADE = T-SUNLIGHT;

ACCESS = elevs_cube > elevcrit;

ACCESSmean = sum(ACCESS,2)*dt/t(end) * T;


BITPERIMAGE = 7920*6004*12;
BITTOTRANSFER = ACCESSmean * 3e08;

IMAGES = floor(BITTOTRANSFER/BITPERIMAGE);

res = 1.05;
timeperimage = res * 7920 * (R0 + alt)/R0 / v / 1e03;

EXPOSURETIME0 = IMAGES * timeperimage;
EXPOSURETIME1 = 5400 * (T/86400);

EXPOSURETIME = min([EXPOSURETIME0,EXPOSURETIME1]);

%% ENERGY - TYPE/MODE

% TYPE 1 MODE 1
PANELS_VEC_11 = [1 0 0
    -1 0 0
    0 1 0
    0 -1 0
    0 0 1
    0 0 -1];
PANELS_POW_11 = [240; 240; 200; 0; 240; 240];

POWPANELS_11 = zeros(6,Nt);
POWERS_11 = zeros(1,Nt);

for i = 1:Nt
    vec = math.xyz2tsw(rsun(:,i),eph(i,4),eph(i,5),eph(i,6))';
    vec = repmat(vec,[6,1]);
    cosfi = dot(PANELS_VEC_11,vec,2);
    cosfi = cosfi.*heaviside(cosfi);
    POWPANELS = PANELS_POW_11.*cosfi;
    POWPANELS_11(:,i) = POWPANELS;
    POWERS_11(i) = sum(POWPANELS);
    if ~lit(i)
        POWPANELS_11(:,i) = 0;
        POWERS_11(i) = 0;
    end
end


plot(t,POWPANELS_11); hold on; grid;
plot(t,POWERS_11);
legend('+T','-T','+S','-S','+W','-W','SUM');



% for i = 1:Nt
%     vec = 2*rand(3,1)-1; vec = vec'/norm(vec);
%     vec = repmat(vec,[6,1]);
%     cosfi = dot(PANELS_VEC_11,vec,2);
%     cosfi = cosfi.*heaviside(cosfi);
%     POWPANELS = PANELS_POW_11.*cosfi;
%     POWERS_11(i) = sum(POWPANELS);
% end

MEANP_11 = mean(POWERS_11(POWERS_11(:) > 0));
MEANP_11_eff = MEANP_11 * 0.85 * 0.99^lifetime;




