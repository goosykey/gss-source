%% TEST SCRIPT FOR SIRIUS

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

%% ORBIT PROPAGATION

[ eph ] = astro.j2circ( INI_STATE, t );
eph = permute(eph,[2,3,1]);

[r,v] = math.randv(eph(:,1), eph(:,2), eph(:,5), eph(:,4), eph(:,3), eph(:,6)); % Keplerian -> Cartesian
r = r'; v = v';

%% SUN

rsun = astro.sun(jt)'; rsun1 = rsun;
rsun = rsun ./ repmat(sqrt(sum(rsun.^2,2)),[1,3]);

d = testsight(r,rsun);

%% PANELS IN TSW

q = 1/sqrt(2);

p1 = [q 0 q];  p1 = repmat(p1, [numel(t),1]);
p2 = [-q 0 q]; p2 = repmat(p2, [numel(t),1]);
p3 = [q 0 -q];  p3 = repmat(p3, [numel(t),1]);
p4 = [-q 0 -q]; p4 = repmat(p4, [numel(t),1]);
p5 = [0 1 0]; p5 = repmat(p5, [numel(t),1]);

%% SUN IN TSW

rsuntsw = rsun*0;

for i = 1:numel(t)
    rsuntsw(i,:) = math.xyz2tsw(rsun(i,:)',eph(i,4),eph(i,5),eph(i,6))';
end

%% POWERS

PC = 2.3; %W

pow1 = PC .* dot(p1, rsuntsw, 2) .* heaviside(dot(p1, rsuntsw, 2)) .* d;

pow2 = PC .* dot(p2, rsuntsw, 2) .* heaviside(dot(p2, rsuntsw, 2)) .* d;

pow3 = PC .* dot(p3, rsuntsw, 2) .* heaviside(dot(p3, rsuntsw, 2)) .* d;

pow4 = PC .* dot(p4, rsuntsw, 2) .* heaviside(dot(p4, rsuntsw, 2)) .* d;

pow5 = PC/3 .* dot(p5, rsuntsw, 2) .* heaviside(dot(p5, rsuntsw, 2)) .* d;

POW = pow1+pow2+pow3+pow4+pow5;

%% PLOTTING

plot(t, pow1); hold on; grid;
plot(t, pow2);
plot(t, pow3);
plot(t, pow4);
plot(t, pow5);

plot(t, POW, 'r', 'linewidth', 2);

legend('Pan1','Pan2','Pan3','Pan4','Pan5','Psum');

%% GROUND

GS = [45 39];
[gscart] = astro.gs_cube(GS, t, jd0);
gscart = permute(gscart,[2 3 1]);


%% CONSUMPTION

P_com = access2(r, gscart, 15) * 5;

P_cam = access2(r, gscart, 61) * 1 .* d;

P_obc  = t'*0 + 0.5; % obc + gps
P_adcs = t'*0 + 0.75;


W_ac_0 = 37; % Wh
W_ac = t*0 + W_ac_0;

charge_eff = 0.95;

for i = 2:numel(t)
    CONS = dt*(P_obc(i) + P_adcs(i) + P_com(i) + P_cam(i)) / 3600;
    CONS = CONS * 1.1; % power supply systemm consumption
    
    W_ac(i) = W_ac(i-1) - CONS / charge_eff + dt*POW(i)/3600 * charge_eff;
    
    if W_ac(i) > W_ac_0
        W_ac(i) = W_ac_0;
    end
end

figure; 
plot(t,W_ac./W_ac_0*100,'r','linewidth',2); hold on; grid;





