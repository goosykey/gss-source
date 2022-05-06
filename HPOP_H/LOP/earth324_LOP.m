function Dy = earth324_LOP(t, y, varargin)
%>>>MAIN MOTION MODEL <<<
%   defines perturbations and motion equations
%
%
%     INPUTS:
%     t : nominal orbit parameters
%     y : single ephemeris point
%       
% 
%     OUTPUTS:
%     Dy : ephemeris DERIVATIVES
% 


%% VARARGIN DEFAULTS

jd0 = astro.jday(2020, 3, 20, 12, 0, 0);

C_ba = 0.0082; %m2/kg
Cxa = 2.2;

F107 = 150;
F81 = 150; 
Kp = 3;


%% PARSE VARARGIN

while ~isempty(varargin)
    switch lower(varargin{1})
        case 'jd0'
            jd0 = varargin{2};
        case 'c_ba'
            C_ba = varargin{2};
        case 'f107'
            F107 = varargin{2};
        case 'f81'
            F81 = varargin{2};
        case 'kp'
            Kp = varargin{2};
        case 'cxa'
            Cxa = varargin{2};
        otherwise
            %varargin{1}
            error(['Unexpected option: ' varargin{1}])
    end
    varargin(1:2) = [];
end

%% INITIALIZE

Dy = y(:);

rxyz = [y(1), y(2), y(3)]';
vxyz = [y(4), y(5), y(6)]';

r = norm(rxyz);
v = norm(vxyz);

Ixyz = cross(rxyz,vxyz); % orbital momentum vector
I = norm(Ixyz);

Om = atan2(Ixyz(1),-Ixyz(2));
in = acos(Ixyz(3)/I);
u = atan2(rxyz(3)/sin(in), rxyz(1)*cos(Om)+rxyz(2)*sin(Om));

jd = jd0 + t/86400;
[month, day, ~] = astro.gdate (jd);
DOY = 30*(month-1) + round(day);

%fprintf('t = %3.2f days ; h = %3.2f km \n',t/86400,r-6378.14);


%% CONSTANTS
del = 66.07354399E03; % simplified J2 (all constants packed)
mu = 398600.4415; % Earth grav.param. km3/s2


%% SUNMOON

[~,rtascsun,declsun] = astro.sun ( jd ); %Vallado



%% nonspherics

recef = math.R3(astro.gstime(jd))*rxyz;

Tj2 = -mu*del/r^4 * sin(2*u) * (sin(in))^2;
Sj2 = mu*del/r^4 * (3*(sin(u))^2*(sin(in))^2-1);
Wj2 = mu*del/r^4 * sin(u) * sin(2*in);
fgxyz = math.tsw2xyz([Tj2;Sj2;Wj2],Om,in,u);

%fgxyz = [0;0;0];


%% ATMOSPHERIC DRAG

atm = 'jr';

if strcmpi(atm, 'gost')
    rho = astro.atmgost2(recef, r-6378.14,DOY,astro.gstime(jd),rtascsun,declsun,F107,F81,Kp); % GOST
elseif strcmpi(atm, 'j77')
    rho = astro.atmJ77(jd,rxyz,F107,F81,Kp); % Jacchia 1977
elseif strcmpi(atm, 'jr')
    rho = astro.atmJR(jd,rxyz,F107,F81,Kp); % Jacchia-Roberts
else % just in case
    rho = astro.atmgost2(recef, r-6378.14,DOY,astro.gstime(jd),rtascsun,declsun,F107,F81,Kp); % GOST
end

Tatm = -Cxa*rho*(v*1000)^2*C_ba/2 * 10^(-3);
fatm = math.tsw2xyz([Tatm; 0; 0],Om,in,u);


%% TOTAL PERTURBATIONS

fx = fgxyz(1) + fatm(1);
fy = fgxyz(2) + fatm(2);
fz = fgxyz(3) + fatm(3);


%% FINALLY -- EQUATIONS

Dy(1) = vxyz(1);
Dy(2) = vxyz(2);
Dy(3) = vxyz(3);

Dy(4) = -mu/r^3*rxyz(1) + fx;
Dy(5) = -mu/r^3*rxyz(2) + fy;
Dy(6) = -mu/r^3*rxyz(3) + fz;

Dy(7) = 0;



end