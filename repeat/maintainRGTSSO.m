%% delta v for maintaining Repeat Groundtrack
% ref. Vallado p. 869


%% CONSTANTS

R0 = 6378.14; %km
mu = 398600.44; %km^3/s^2
deg = pi/180;
J2 = 1.0826267E-03;
omega = 7.2921151467e-5;     % Earth inertial rotation rate (rad/sec)

%% INITIAL

ALT = 560.993907; %km

corridor = 0.1*deg;

a = R0 + ALT;
e = 0;
i = 97.640225*deg;
Om = 72.06*deg;

STARTTIME = [2018 1 1 0 0 0];

jd0 = astro.jday(STARTTIME(1), STARTTIME(2), STARTTIME(3), STARTTIME(4), STARTTIME(5), STARTTIME(6));

%% IMPLEMENTATION

p = a * (1 - e^2); % focal parameter
n = sqrt(mu/a^3); % keplerian mean motion
RAANdot = -3*n*R0^2*J2/2/p^2 * cos(i);

dPda = 3*pi/n/a * (1 + .5*J2*(R0/a)^2*(4*(cos(i))^2-1));
dPdi = 12*pi/n * J2 * (R0/a)^2 * sin(2*i);
dRAANdotda = -3.5 * RAANdot / a;
dRAANdotdi = -RAANdot * tan(i);

Pkepl = 2*pi * a * sqrt(a / mu); % keplerian period
P = 2*pi/n * (1 - 3*J2/2 * (R0/a)^2 * (3 - 4*sin(i)^2)); % nodal period

dlarevnom = (omega - RAANdot) * P;
dlarevda = R0 * (omega - RAANdot) * dPda - dRAANdotda * R0 * P;
dlarevdi = R0 * (omega - RAANdot) * dPdi - dRAANdotdi * R0 * P;

% y = [a 0.0005 0.0005 Om i -3*pi/4 200];
% Dy = earth353_FULL(jd0, y);
% 
% dadt = Dy(1); % only for small e !!!
% didt = Dy(5);

dadt = -6.2063e-08; %km/sec
didt = -1.515E-10; %rad/sec

k2 = 1/P * (dlarevda*dadt + dlarevdi*didt);
k1 = sqrt(2*k2 * (-corridor));

tdrift = -k1/k2;

deltaa = k1 * P / dlarevda;

deltav = n/2 * deltaa;
fprintf('\n    <here we go>  \n\n');

fprintf('Period Keplerian =     %3.2f sec \n',Pkepl);
fprintf('Period odal     =     %3.2f sec \n\n',P);

fprintf('dlanom =     %12.6f deg \n',dlarevnom/deg);
fprintf('Tdrift =     %12.2f sec \n',tdrift);
fprintf('Delta a =    %12.6f m \n',deltaa * 1e3);
fprintf('Delta V =    %12.6f m/s \n',deltav * 1e3);