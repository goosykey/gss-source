function m = GetMagnitude(varargin)
%% References

% The function returns magtidue of the Satellite's reflector

% 1. Canady Jr, John E., and John L. Allen Jr.
% "Illumination from space with orbiting solar-reflector spacecraft." (1982).
% 2. Allen, Clabon Walter. "Astrophysical quantities." (1973).
% 3. Hottel, Hoyt C. "A simple model for estimating the transmittance of direct solar radiation through clear atmospheres." Solar energy 18.2 (1976): 129-134.

%% PRE-DEFINE DEFAULT VALUES

ro = 0.92; % mirror reflectance Maylar or Kapton coated with aluminum [1,2]
solarflux = 1360;  % Solar flux intensity, W/m2
sm = 80; % mirror area, m^2
alpha = 5.21; % angle of scattering from mirror surface, degrees
gamma = 120; % is the angle between incident and reflected light beam, degrees
elev = 80; % elevation angle, degrees
elev(elev(:)<0) = nan;
sat2spotdist = 550e3; % distance from the mirror to the spot, m

%% VARARGIN CHECK

nvarargs = length(varargin);

if rem(nvarargs,2) ~= 0
    error('Invalid number of arguments' )
end

while ~isempty(varargin)
    switch lower(varargin{1})
        case 'ro'
            ro = varargin{2};
        case 'solarflux'
            solarflux = varargin{2};
        case 'sm'
            sm = varargin{2};
        case 'gamma'
            gamma = varargin{2};
        case 'elev'
            elev = varargin{2};
            elev(elev<=0) = NaN;
        case 'alpha'
            alpha = varargin{2};
        case 'sat2spotdist'
            sat2spotdist = varargin{2};
        otherwise
            error(['Unexpected option: ' varargin{1}])
    end
    varargin(1:2) = [];
end

%% Calculations

tau = 0.1283 + 0.7559*exp(-0.3878*secd(90-elev)); % atmospheric transmissivity [3]

I = (solarflux.*ro.*tau.*sm.*cosd(gamma/2).*sind(elev))./(pi*(sat2spotdist.*tand(alpha/2)).^2);

m = -2.5*log10(I*683./2.54e-6);

end