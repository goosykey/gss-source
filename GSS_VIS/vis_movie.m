function fig = vis_movie(t, ephcart, varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%   ephcart EXPECTED CUBE (Nsat, Nt, 6)

%% CONST
deg = pi/180;
erot    = 7.2921158553e-5; % earth rotation rate (radians/sec)

timeres = 0.02; % THIS IS REAL TIME RESOLUTION

%% PRE-DEF

timescale = 240; % X REAL TIME
jd0 = 2451545;
satstyles = 'or';
satcolors = [0.75 0.75 0.75];

orbitpersistence = 43200;

%% VARARGIN CHECK

while ~isempty(varargin)
    switch lower(varargin{1})
        case 'jd0'
            jd0 = varargin{2};
        case 'satstyles'
            satstyles = varargin{2};
        case 'satcolors'
            satcolors = varargin{2};
        case 'timescale'
            timescale = varargin{2};
        otherwise
            error(['Unexpected option: ' varargin{1}])
    end
    varargin(1:2) = [];
end

%%DEFINITIONS & CALCULATIONS

[Nsat, ~, ~] = size(ephcart);

dt = t(2) - t(1);

dtnew = timeres * timescale; % model time step while 0.04 real seconds pass
tnew = t(1):dtnew:t(end);
Ntnew = numel(tnew);

ephcartnew = zeros(Nsat,numel(tnew),3);

for i = 1:3
    ephcartnew(:,:,i) = interp1(t,ephcart(:,:,i)',tnew)'; % NEW EPH DATA
end

globerotangle = erot*dtnew/deg; % IN DEGREES

%% IMPL

[ fig, globe, misc ] = vis_earthdraw(jd0, 'quality','low');

X = ephcartnew(:,1,1); Y = ephcartnew(:,1,2); Z = ephcartnew(:,1,3);
sats = plot3(X,Y,Z,satstyles, 'markerfacecolor', satcolors);

perssteps = round(orbitpersistence/dtnew);

X = ephcartnew(:,1:2,1); Y = ephcartnew(:,1:2,2); Z = ephcartnew(:,1:2,3);
orbits = plot3(X',Y',Z','-r');

for i = 2:Ntnew
    rotate([globe; misc],[0 0 1], globerotangle);
    sats.XData = ephcartnew(:,i,1);
    sats.YData = ephcartnew(:,i,2);
    sats.ZData = ephcartnew(:,i,3);
    
    if mod(i,1000) == 0 || i == 2 % Plot orbits
        if i + perssteps/2 <= Ntnew
            if i - perssteps/2 >= 1
                persrange = round(i - perssteps/2) : round (i + perssteps/2);
            else
                persrange = 1 : round (i + perssteps/2);
            end
        else
            if i - perssteps/2 >= 1
                persrange = round(i - perssteps/2) : Ntnew;
            else
                persrange = 1 : Ntnew;
            end
        end
    end
    X = ephcartnew(:,persrange,1); Y = ephcartnew(:,persrange,2); Z = ephcartnew(:,persrange,3);
    delete(orbits);
    orbits = plot3(X',Y',Z','-r');
    pause (timeres);
end


end

