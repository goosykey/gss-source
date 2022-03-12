function [ masses, times ] = iophboost( N, m_disp, plotornot )
%IOPHASING_BOOST phasing by rocket booster
%   Phasing in-orbit by Rocket-Booster, for NF @ TECHCOM 11.05.17

if nargin < 2
    m_disp = 500;
    plotornot = 1;
elseif nargin < 3
    plotornot = 1;
end

Nsc = 14;
m_sc = 300;
alt = 650;
R0 = 6378.14;
mu = 398600.44;
deg = pi/180;
Isp = 333.2*9.8; % m/s

du = 360*deg/Nsc; % radians

m_dry = 980; % dry mass fregat, kg
m_debris = m_dry + m_disp + 100;

r = R0+alt;
T1 = 2*pi*sqrt(r^3/mu);

t = T1/2/pi*du; % seconds

T2 = T1 + t/N;
a2 = (sqrt(mu)*T2/2/pi)^(2/3);
rp = r;
ra = 2*a2-rp;

h2 = sqrt(2*mu) * sqrt(ra*rp/(ra+rp));
h1 = sqrt(mu*r);

delta_v = h2/r - h1/r; % km/s
time = T2*N;

m_after = zeros(Nsc,1);
m_before = zeros(Nsc,1);

M2 = m_debris;

for i = Nsc:-1:1
    m_after(i) = M2; % after jettison
    M2 = M2 + m_sc; % jettison mass
    m_before(i) = M2; % before jettison
    M1 = M2 * exp(2*delta_v*1e3/Isp); % after jettison
    M2 = M1;
end

times = zeros(2*Nsc,1);
masses = zeros(2*Nsc,1);

for i = 0:2:(2*Nsc-2)
    times(i+1) = i/2*time;
    times(i+2) = i/2*time;
    masses(i+1) = m_before(i/2+1);
    masses(i+2) = m_after(i/2+1);
end

if ~plotornot
    return
end

s = sprintf('N = %g',N);
plot(times/3600,masses,'o-', 'Displayname',s);
hold on;
grid on;
xlabel('time, hours');
ylabel('cluster mass, kg');
legend('show');


end

