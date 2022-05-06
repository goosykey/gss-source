function [rho] = atmJR(jd, Rxyz, F107, F81, kp)
%ATMJR Summary of this function goes here
%   Vallado p. 1001

% clear;clc;
% jd = 2458599.125;
% Rxyz = [4.7774e+03 ; -295.9764; 4.2154e+03] * ((6378.14+600)/6378.14);
% F107 = 150;
% F81 = 150;
% kp=3;


%% CONST

R0 = 6378.137;
deg = pi/180;

f = 1/298.257223563; % ellipsoid flattening
e2 = 2*f-f^2; % ellipsoid first eccentricity squared
Rp = (1-f)*R0;

R = 8.31432; % J/mole*K, universal gas constant
g_SL = 9.80665; % free-fall acceleration @ sea level

Av = 6.02257e23; % Avogadro

eps = 23.4*deg; % obliquity of the ecliptic



%% LOCAL COORDS

rI = Rxyz(1);
rJ = Rxyz(2);
rK = Rxyz(3);

%% SUN

[rsun, ~, declsun] = astro.sun(jd);

rx = rsun(1);
ry = rsun(2);


%% EVALUATING TEMPERATURE

% Nighttime global exospheric temperature
T_c = 379 + 3.24 * F81 + 1.3 * (F107-F81);

% coefs from Vallado p. 1001
LHA_sun = (rx*rJ - ry*rI) / abs(rx*rJ - ry*rI) ...
    * acos( (rx*rI + ry*rJ)/sqrt(rx^2+ry^2)/sqrt(rI^2+rJ^2) );
fi_gd = atan(1/(1-f)^2 * rK / sqrt(rI^2 + rJ^2));
eta = abs(fi_gd - declsun) / 2;
theta = abs(fi_gd + declsun) / 2;
tau = LHA_sun - 37*deg + 6*deg * sin(LHA_sun + 43*deg);

%% ELLIPSOIDAL PARAMETERS

% auxiliary quantities (Vallado p.138 3-7)
C0 = R0/sqrt(1-e2*sin(fi_gd)^2);
S0 = R0*(1-e2)/sqrt(1-e2*sin(fi_gd)^2);
b = (C0*cos(fi_gd)^2 + S0*sin(fi_gd)^2);
c = C0*C0*cos(fi_gd)^2 + S0*S0*sin(fi_gd)^2 - norm(Rxyz)^2;
h_ellp = -b + sqrt(b^2 - c);

%% EVALUATING TEMPERATURE (continues)

% uncorrected Exospheric temperature
T_unc = T_c * (1 + 0.3 * (sin(theta)^2.2 + (cos(eta)^2.2-sin(theta)^2.2)*cos(tau/2)^3));

deg1 = 1;

if h_ellp > 200
    dT_corr = 28*deg1 * kp + 0.03*exp(kp);
else
    dT_corr = 14*deg1 * kp + 0.02*exp(kp);
end

% corrected exospheric temperature
T_corr = T_unc + dT_corr;

% inflection point temperature
Tx = 371.6678*deg1 + 0.0518806 * T_corr - 294.3505*deg1 * exp(-0.00216222*T_corr);

T0 = 183;

C = [-89284375.0 3542400.0 -52687.5 340.5 -0.8];

%T90 required later for Roberts' corrections to density
T90 = Tx;
for i = 0:4
        T90 = T90 + (Tx-T0)/35^4 * C(i+1) * 90^i;
end

if h_ellp < 125
    T = Tx;
    for i = 0:4
        T = T + (Tx-T0)/35^4 * C(i+1) * h_ellp^i;
    end
else
    T = Tx + 2/pi * (T_corr-Tx) * atan(0.95*pi*((Tx-T0)/(T_corr-Tx))*((h_ellp-125)/35)*(1+4.5e-06*(h_ellp-125)^2.5));
    
end


%% ROBERT'S CORRECTIONS TO TEMPERATURE

l_Jacchia = 12315.3554; % jacchia

l_i_Draper = [0.1031445e05 0.2341230e01 0.1579202e-02 -0.1252487e-05 0.2462708e-09];

l_Draper = 0;

for i = 0:4
    l_Draper = l_Draper + l_i_Draper(i+1) * T_corr^i;
end

% in the following correction, use either l_Jacchia (by Jacchia) or l_Draper
l = l_Draper;

if h_ellp >= 125
    T = T_corr - (T_corr-Tx) * exp(-((Tx-T0)/(T_corr-Tx))*((h_ellp-125)/35)*(l/(Rp+h_ellp)));
end

%% EVALUATING DENSITY

% Vallado p.567 table 8-4
h_base = [0 25 30 40 50 60 70 80 90 100 110 120 130 140 150 180 200 250 300 350 400 450 500 600 700 800 900 1000 inf];
rho_base = [1.225 3.899e-02 1.774e-02 3.972e-03 1.057e-03 3.206e-04 ...
    8.770e-05 1.905e-05 3.396e-06 5.297e-07 9.661e-08 2.438e-08 8.484e-09 ...
    3.845e-09 2.070e-09 5.464e-10 2.789e-10 7.248e-11 2.418e-11 9.518e-12 ...
    3.725e-12 1.585e-12 6.967e-13 1.454e-13 3.614e-14 1.170e-14 5.245e-15 ...
    3.019e-15];
H_base = [7.249 6.349 6.682 7.554 8.382 7.714 6.549 5.799 5.382 5.877 7.263 9.473 12.636 16.149 ...
    22.523 29.740 37.105 45.546 53.628 53.298 58.515 60.828 63.822 71.835 88.667 124.64 181.05 268.00];

% index of current h_ellp in table 8-4
[~,index] = max(diff(h_ellp < h_base));

rho_0 = rho_base(index);
h_0 = h_base(index);
H = H_base(index);

% standard exponential density, Vallado 8-33 and B-6
rho_std = rho_0 * exp(-(h_ellp-h_0)/H);

% geomagnetic effect on density
if h_ellp < 200 
    delta_log10_rho_G = 0.012 * kp + 1.2e-05 * exp(kp);
else
    delta_log10_rho_G = 0;
end

JD1958 = jd - 2436204;

T1958 = JD1958/365.2422;

% seasonal latitudinal variation
delta_log10_rho_LT = 0.014*(h_ellp-90) * sin(2*pi*T1958+1.72)  ...
    * sin(fi_gd) * abs(sin(fi_gd)) * exp(-0.0013*(h_ellp-90)^2);

tau_SA = T1958 + 0.09544 * ((0.5 + 0.5 * sin(2*pi*T1958 + 6.035))^1.65 - 0.5);

% semi-annual variations
delta_log10_rho_SA = (5.876e-07*h_ellp^2.331+0.06328) * exp(-0.002868*h_ellp) ...
    * (0.02835 + (0.3817 + 0.17829 * sin(2*pi*tau_SA+4.137)) ...
    * sin(4*pi*tau_SA+4.259));

% total correction
delta_log10_rho_corr = delta_log10_rho_G + delta_log10_rho_LT + delta_log10_rho_SA;

% final corrected density
rho = rho_std * 10^(delta_log10_rho_corr); % Jaccha 1971 density

%% ROBERT'S CORRECTION TO DENSITY

Cstar = C*0;

Cstar(1) = 35^4*Tx / (C(5)*(Tx-T0));
for i = 2:5
    Cstar(i) = C(i)/C(5);
end


roots_all = roots(Cstar(end:-1:1));
roots_real = roots_all(imag(roots_all)==0);
root_pimag = roots_all(imag(roots_all) >0);

xr1 = real(roots_real(1));
xr2 = real(roots_real(2));
xr3 = real(root_pimag);
xi3 = imag(root_pimag);

Xstar = -2*xr1*xr2*Rp*(Rp^2 + 2*xr3*Rp + xr3^2 + xi3^2);
V = (Rp+xr1)*(Rp+xr2)*(Rp^2 + 2*xr3*Rp + xr3^2 + xi3^2);
U_xr1 = (xr1+Rp)^2*(xr1^2 - 2*xr3*xr1 + xr3^2 + xi3^2) * (xr1-xr2);
U_xr2 = (xr2+Rp)^2*(xr2^2 - 2*xr3*xr2 + xr3^2 + xi3^2) * (xr1-xr2);
W_xr1 = xr1*xr2*Rp*(Rp+xr1)*(Rp+(xr3^2 + xi3^2)/xr1);
W_xr2 = xr1*xr2*Rp*(Rp+xr2)*(Rp+(xr3^2 + xi3^2)/xr2);

alpha = [3144902516.672729 -123774885.4832917 1816141.096520398 ...
    -11403.31079489267 24.36498612105595 0.008957502869707995];
beta = [-52864482.17910969 -16632.50847336828 -1.308252378125 0 0 0];

B = alpha + beta * Tx/(Tx-T0);

S_xr1 = 0; S_xr2 = 0; S_minusRp = 0;

for i = 0:5
    S_xr1 = S_xr1 + B(i+1)*xr1^i;
    S_xr2 = S_xr2 + B(i+1)*xr2^i;
    S_minusRp = S_minusRp + B(i+1)*(-Rp)^i;
end

p2 = S_xr1/U_xr1;
p3 = -S_xr2/U_xr2;
p5 = S_minusRp/V;
p4 = 1/Xstar * ( B(1) - xr1*xr2*Rp^2 * (B(5)+B(6)*(2*xr3+xr1+xr2-Rp)) ...
    + W_xr1*p2 - xr1*xr2*B(6)*Rp*(xr3^2 + xi3^2) + W_xr2*p3 ...
    + xr1*xr2*(Rp^2-xr3^2-xi3^2)*p5);
p6 = B(5) + B(6)*(2*xr3+xr1+xr2-Rp) - p5 - 2*(xr3+Rp)*p4 ...
    - (xr2+Rp)*p3 - (xr1+Rp)*p2;
p1 = B(6) - 2*p4 - p3 - p2;

A = [-435093.363387 28275.5646391 -765.33466108 11.043387545 ...
    -0.08958790995 0.00038737586 0.000000697444];

f1 = 35^4*Rp^2/C(5);

F1 = ((h_ellp+Rp)/(90+Rp))^p1 * ((h_ellp-xr1)/(90-xr1))^p2 ...
    * ((h_ellp-xr2)/(90-xr2))^p3 ...
    * ((h_ellp^2-2*xr3*h_ellp+xr3^2+xi3^2)/(8100-180*xr3+xr3^2+xi3^2))^p4;
F2 = (h_ellp-90)*(f1*A(7) + p5/(h_ellp+Rp)/(90+Rp)) ...
    + p6/xi3 * atan(xi3*(h_ellp-90)/(xi3^2+(h_ellp-xr3)*(90-xr3)));



M_h_ellp = 0;
for i = 0:6
    M_h_ellp = M_h_ellp + A(i+1) * h_ellp^i;
end

k = -g_SL/ R / (Tx-T0);

rho90 = 3.46e-06;
M90 = 28.82678;

Ms = 28.96;
M = [28.0134 39.948 4.0026 31.9988 15.9994 1.00797];
a = [0 0 -0.38 0 0 0];
mu = [0.78110 0.0093432 0.61471e-05 0.161778 0.095544 0];

if h_ellp < 150 % TEMP
    h_ellp = 150.01;
end

if h_ellp > 90
    if h_ellp <= 100
        T_Roberts_125plus = T_corr - (T_corr-Tx) * exp(-((Tx-T0)/(T_corr-Tx))*((h_ellp-125)/35)*(l/(Rp+h_ellp)));
        rho_std = rho90*T90/M90 * M_h_ellp/T_Roberts_125plus * F1^k * exp(k*F2);
    elseif h_ellp <= 125 % in 100..125 we always perform calculations for rho or rho125

        q2 = 1/U_xr1;
        q3 = -1/U_xr2;
        q5 = 1/V;
        q4 = (1 + xr1*xr2*(Rp^2-xr3^2-xi3^2)*q5 + W_xr1*q2 + W_xr2*q3)/Xstar;
        q6 = -q5 - 2*(xr3+Rp)*q4 - (xr2+Rp)*q3 - (xr1+Rp)*q2;
        q1 = -2*q4-q3-q2;
        
        F3 = ((h_ellp+Rp)/(100+Rp))^q1 * ((h_ellp-xr1)/(100-xr1))^q2 ...
            * ((h_ellp-xr2)/(100-xr2))^q3 ...
            * ((h_ellp^2-2*xr3*h_ellp+xr3^2+xi3^2)/(10000-200*xr3+xr3^2+xi3^2))^q4;
        
        F4 = (q5*(h_ellp-100)/(h_ellp+Rp)/(Rp+100)) ...
            + q6/xi3*atan(xi3*(h_ellp-100)/(xi3^2+(h_ellp-xr3)*(100-xr3)));
        
        T100 = Tx - 0.94585589 * (Tx - T0);
        
        xsi = [0.1985549e-10 -0.183349e-14 0.1711735e-17 -0.1021474e-20 ...
            0.3727894e-24 -0.7734110e-28 0.7026942e-32];
        
        rho100 = 0;
        for i = 0:6
            rho100 = rho100 + Ms * xsi(i+1) * T_corr^i;
        end
        
        
        
        rho_std = 0;
        for i = 1:5
            rho_std = rho_std + rho100*M(i)/Ms*mu(i) ...
                * (T100/T)^(1+a(i)) ...
                * F3^M(i)*k*f1 * exp(M(i)*k*f1*F4);
        end
        
    else % above 125 km (THIS IS THE ONLY PART THAT ACTUALLY WORKS BUT WE DON'T REALLY NEED THE REST)
        % TABLE B-4 Vallado p.1010
        delta_ij = [0.1093155e02 0.8049405e01 0.7646886e01 0.9924237e01 0.1097083e02
            0.1186783e-02 0.2382822e-02 -0.4383486e-03 0.1600311e-02 0.6118742e-04
            -0.1677341e-05 -0.3391366e-05 0.4694319e-06 -0.2274761e-05 -0.1165003e-06
            0.1420228e-08 0.2909714e-08 -0.2894886e-09 0.1938454e-08 0.9239354e-10
            -0.7139785e-12 -0.1481702e-11 0.9451989e-13 -0.9782183e-12 -0.3490739e-13
            0.1969715e-15 0.4127600e-15 -0.1270838e-16 0.2698450e-15 0.5116298e-17
            -0.2296182e-19 -0.4837461e-19 0.0 -0.3131808e-19 0.0];
        
        rho125_i = [0 0 0 0 0 0];
        for i = 1:5
            summ = 0;
            for j = 0:6
                summ = summ + delta_ij(j+1,i) * T_corr^j;
            end
            rho125_i(i) = M(i)*10^summ/Av * 1e3;
            % in the original expression for rho125 in vallado Avogadro is
            % missed, besides, it seems to be in g/cc
        end
        
        gamma = M*g_SL*Rp^2/R/l/T_corr * ((T_corr-Tx)/(Tx-T0))*35/6481.766;
        rho_i = rho125_i .* (Tx/T).^(1+a+gamma) ...
            .* ((T_corr-T)/(T_corr-Tx)).^gamma;
        
        % seasonal variation of Helium
        delta_log10_rho_He = 0.65 * abs(declsun/eps) * (sin(pi/4-fi_gd*declsun/2/abs(declsun))^3-0.35355);
        rho_i(3) = rho_i(3) * 10^delta_log10_rho_He;
        
        % account for Hydrogen
        T500 = T_corr - (T_corr-Tx) * exp(-((Tx-T0)/(T_corr-Tx))*((500-125)/35)*(l/(Rp+500)));
        rho500_H = M(6)/Av * 10^(73.13-(39.4-5.5*log10(T500))*log10(T500));
        if h_ellp >= 500
            rho_i(6) = rho500_H*(T500/T)^(1+a(6)+gamma(6)) ...
                * ((T_corr-T)/(T_corr-T500))^gamma(6);
        else
            rho_i(6) = 0;
        end
        
        rho_std = sum(rho_i);
    end
end

rho = rho_std;

end