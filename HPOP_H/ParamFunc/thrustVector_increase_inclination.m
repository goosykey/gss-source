function [thr] = thrustVector_increase_inclination(jd, Y)
%thrustVector_increase_inclination THRUST VECTOR FOR INCREASE INCLINATION
%   This function returns thrust direction in the TSW frame, always along W
%   or -W, always corresponding to increasing inclination.
% 
%   INPUT:
%       jd : epoch in Julian Day (non-modified)
%       y  : phase vector in earth353 format:
%          [p
%           lambda1
%           lambda2
%           Omega
%           inc
%           u]
%
%   OUTPUT:
%       thr : thrust direction in TSW


%% INITIALIZE



%% CONSTANTS

deg = pi/180;


%% IMPLEMENTATION

[~,~,~,~,~,u] = math.cart2kep(Y(1), Y(2), Y(3), Y(4), Y(5), Y(6));

u = mod (u, 2*pi);

if u <= 90*deg || u > 270*deg
    thr = [0;0;-1];
else
    thr = [0;0;1];
end


end

