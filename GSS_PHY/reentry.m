function [ tact,yact,tfall,yfall ] = reentry( jd0, y0 )
%REENTRY Summary of this function goes here
%   Detailed explanation goes here

global SATDATA;

SATDATA{3} = ones(360,1) * pi;
flg = zeros(360,1); flg(1:45) = 1; flg(316:360) = 1; flg(:)=1;
SATDATA{4} = boolean(flg);

y0(7) = y0(7) - 3;

h = 100;
options = odeset('MaxStep', 240, 'Events',@eventa);

[tact, yact] = ode45(@earth353_FULL, (jd0*86400):h:(jd0*86400 + 6e7), y0, options);

fprintf('Active stage fin\n');

% SATDATA{4} = zeros(360,1);
% options = odeset('MaxStep', 240, 'Events',@eventb);
% 
% [tfall, yfall] = ode45(@earth353_FULL, (tact(end)):h:(tact(end) + 6e7), yact(end,:), options);

tfall = 0; yfall = 0;

end

function [value,isterminal,direction] = eventa( t, y )
%ODE_EV Summary of this function goes here
%   Detailed explanation goes here

%fprintf('Si\n');

Y = math.unorbitize(y');
[R,~] = math.randv(Y(1),Y(2),Y(5),Y(4),Y(3),Y(6)-Y(3));

R = sqrt(sum(R.^2));

value = R - 6378.14 - 200;     % detect height = 0
isterminal = 1;   % stop the integration
direction = 0;   % negative direction

end

function [value,isterminal,direction] = eventb( t, y )
%ODE_EV Summary of this function goes here
%   Detailed explanation goes here


Y = math.unorbitize(y');
[R,~] = math.randv(Y(1),Y(2),Y(5),Y(4),Y(3),Y(6)-Y(3));

R = sqrt(sum(R.^2));

value = R - 6378.14 - 50;     % detect height = 0
isterminal = 1;   % stop the integration
direction = 0;   % negative direction

end