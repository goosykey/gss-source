R0 = 6378.14;
deg = pi/180;
mu0 = 398600.44;

jd0 = astro.jday(2018, 3, 20, 16, 15, 0);

%% Define satellites

t = 0:60:7200;
jt = jd0 + t/86400;
Nt = numel(t);

alt = 600;

a = R0+alt;
e = 0;
om = 0;
in = astro.sunsyn(alt)*deg;

ini0 = [a e om 0 in 0
    a e om 1 in 1
    a e om 2 in 2
    a e om 3 in 3];

Nsat = numel(ini0(:,1));

T = 2*pi*sqrt(ini0(:,1).^3/mu0);


%% Define ground stations

GS = [53 37
    0 0
    40 -135];

Ngs = numel(GS(:,1));


%% Get coordinates

[ eph ] = astro.j2circ( ini0, t );

ephcart = eph .* 0;

for i = 1:Nsat
    [R,V] = math.randv(eph(i,:,1),eph(i,:,2),eph(i,:,5),eph(i,:,4),eph(i,:,3),eph(i,:,6)-eph(i,:,3));
    ephcart(i,:,1:3) = R';
    ephcart(i,:,4:6) = V';
end

gscart = zeros(Ngs,Nt,3);

for i = 1:Ngs
    X = R0 * cosd(GS(i,1)) .* cos(GS(i,2)*deg - astro.gstime(jt));
    Y = R0 * cosd(GS(i,1)) .* sin(GS(i,2)*deg - astro.gstime(jt));
    Z = R0 * sind(GS(i,1)) + jt*0;
    gscart(i,:,:) =[X; Y; Z]';
end

%% BULLSHIT

K = 1:Nsat;

% Nrev
Nrev = ceil(repmat(t,[Nsat,1]) ./ repmat(T,[1,Nt]));
Nrev(Nrev==0) = 1;


%% FIN

clear X Y Z R V








