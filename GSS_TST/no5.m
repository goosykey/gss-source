%% GEOMETRY INPUT

R0 = 6378.14;
deg = pi/180;
alt = 1000;
R = R0+alt;

%% COMMON LINK DATA
L_line = 2; %dB
f = 1.25; %GHz;
d_ss = 2; %m
br = 256e03; %bps
Ts = 10*log10(290); %dBK

losses_introp = 1; % dB/km

Gain_Tx = 17.8 + 20*log10(d_ss) + 20*log10(f); %dB
theta_Tx = 21/f/d_ss; %deg

%% GEOMETRY

th = [51.50, -0.12];
je = [53.41, -2.98];
ss = [47.37, 8.54];

r_th = [R0*cosd(th(1))*cosd(th(2)) R0*cosd(th(1))*sind(th(2)) R0*sind(th(1))];
r_je = [R0*cosd(je(1))*cosd(je(2)) R0*cosd(je(1))*sind(je(2)) R0*sind(je(1))];
r_ss = [R*cosd(ss(1))*cosd(ss(2)) R*cosd(ss(1))*sind(ss(2)) R*sind(ss(1))];

r1_ts = r_ss - r_th;
r1_js = r_ss - r_je;

alpha_th = acosd(dot(r_ss,r1_ts)/norm(r_ss)/norm(r1_ts));
alpha_je = acosd(dot(r_ss,r1_js)/norm(r_ss)/norm(r1_js));

alpha_tj = acosd(dot(r1_ts,r1_js)/norm(r1_ts)/norm(r1_js));

elev_th = acosd(sind(alpha_th)*R/R0);
elev_je = acosd(sind(alpha_je)*R/R0);

dist_th = norm(r1_ts);
dist_je = norm(r1_js);

dist_introp_th = 10 / cosd(alpha_th);
dist_introp_je = 10 / cosd(alpha_je);

%% THERESA LINK

Gain_th = 0;

L_fs_th = 92.45 + 20*log10(dist_th) + 20*log10(f);
L_introp_th = losses_introp * dist_introp_th;

EbN0_th = 9.6; %dB
CN0_th = EbN0_th + 10*log10(br);

GT_th = Gain_th - Ts; %dB
EIRP = CN0_th - GT_th + L_fs_th + L_introp_th + L_line - 228.6;

Ptx = EIRP - Gain_Tx + L_line;
Ptx_W = 10^(Ptx/10);

%% JEREMY LINK

L_line_je = 15;

Gain_je = 19; %dBi

L_fs_je = 92.45 + 20*log10(dist_je) + 20*log10(f);
L_introp_je = losses_introp * dist_introp_je;

GT_je = Gain_je - Ts; %dB

CN0_je = EIRP - 3 + GT_je - L_fs_je - L_introp_je - L_line_je + 228.6;
EbN0_je = CN0_je - 10*log10(br);

D_je = 10^((Gain_je-17.8-20*log10(f))/20);