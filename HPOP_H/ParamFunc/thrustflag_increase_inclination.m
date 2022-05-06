function [thf] = thrustflag_increase_inclination(jd, Y)
%SM_IDEALPANELS CSArea for ideal solar panels
%   This function returns thrust flag to activate thruster burns in the
%   proximity of the orbital nodes to regulate the trade-off point
%   between propellant consumption and manoeuvring efficiency.
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
%       thf : thrust flag (0/1)



%% CONSTANTS

deg = pi/180;

%% INITIALIZE

du = 20*deg; % half-width of the angle corridor near the nodes (max angular distance from a node)




%% IMPLEMENTATION

[~,~,~,~,~,u] = math.cart2kep(Y(1), Y(2), Y(3), Y(4), Y(5), Y(6));

u = mod (u, 2*pi);

if u <= du || u > 2*pi-du || (u >= pi-du && u < pi+du)
    thf = 1;
else
    thf = 0;
end


end