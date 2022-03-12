function [value,isterminal,direction] = ode_ev( L5, y )
%ODE_EV Summary of this function goes here
%   Detailed explanation goes here

global jdsfin;

value = y(6) - jdsfin;     % detect height = 0
isterminal = 1;   % stop the integration
direction = 0;   % negative direction

end

