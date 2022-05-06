classdef OTTSegment < handle
    %OTTSEGMENT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        type, ...
            TState_jd, TState_ifKepl, TState_ini0, ...
            TProp_FM, TProp_PF, TProp_endState, ...
            TiBurn_ifTSW, TiBurn_dV, ...
            TfBurn_ifTSW, TfBurn_dV, ...
            linecolor
    end
    
    methods
        function obj = OTTSegment(varargin)
            obj.type = 'state';
            obj.TState_jd = 2459215.5;
            obj.TState_ifKepl = true;
            obj.TState_ini0 = [6978.137, 0, 0, 0, 1, pi/4, 200];
            obj.TProp_FM = forceModelSet;
            obj.TProp_PF = {};
            obj.TProp_endState = [];
            obj.TiBurn_ifTSW = true;
            obj.TiBurn_dV = [1;0;0];
            obj.linecolor = 'r';
            
            while ~isempty(varargin)
                switch lower(varargin{1})
                    case {'type'}          
                        obj.type = varargin{2};
                    case {'jd'}
                        mustBeNumeric(varargin{2});
                        obj.TState_jd = varargin{2};                        
                    case {'ephtype'}
                        obj.TState_ifKepl = strcmp(varargin{2},'kepler');
                    case {'ini0'}
                        obj.TState_ini0  = varargin{2};
                    case {'forcemodels'}
                        obj.TProp_FM = varargin{2};
                    case {'parameters'}
                        obj.TProp_PF = varargin{2};
                    case {'tsw'}
                        obj.TiBurn_ifTSW = varargin{2};
                    case {'dv'}
                        obj.TiBurn_dV = varargin{2};
                    case {'linecolor'}
                        obj.linecolor = varargin{2};
                    otherwise
                        error(['Unexpected option: ' varargin{1}])
                end
                varargin(1:2) = [];
            end
            
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

 % 
OTT/StopConditions/event_apogee.m
 
function [value, isterminal, direction] = event_apogee(T, Y)
%EVENT_SMA Summary of this function goes here
%   Detailed explanation goes here

@(a)T;

[~,~,om,~,~,u] = math.cart2kep(Y(1), Y(2), Y(3), Y(4), Y(5), Y(6));

nu = u-om;

nu = mod(nu,2*pi);

value = nu-pi;
isterminal = 1;
direction = 1;

end

 % 
OTT/StopConditions/event_asc_node.m
 
function [value, isterminal, direction] = event_asc_node(T, Y)
%EVENT_SMA Summary of this function goes here
%   Detailed explanation goes here

@(a)T;

[~,~,~,~,~,u] = math.cart2kep(Y(1), Y(2), Y(3), Y(4), Y(5), Y(6));

u(u>2*pi) = u-2*pi;

value = u;
isterminal = 1;
direction = 1;

end

 % 
OTT/StopConditions/event_perigee.m
 
function [value, isterminal, direction] = event_perigee(T, Y)
%EVENT_SMA Summary of this function goes here
%   Detailed explanation goes here

@(a)T;

[~,~,om,~,~,u] = math.cart2kep(Y(1), Y(2), Y(3), Y(4), Y(5), Y(6));

nu = u-om;

nu = mod(nu,2*pi);
nu(nu>pi) = nu-2*pi;

value = nu;
isterminal = 1;
direction = 1;

end

 % 
OTT/StopConditions/event_SMA.m
 
function [value, isterminal, direction] = event_SMA(T, Y)
%EVENT_SMA Summary of this function goes here
%   Detailed explanation goes here

@(a)T;

R0 = 6378.14;

p = Y(1);
l1 = Y(2);
l2 = Y(3);
u = Y(6);

e = sqrt(l1.^2 + l2.^2);
om = atan2(l2,l1);
nu = u-om;
r = p / (1 + e*cos(nu));

value = (r < (R0+150));

isterminal = 1;

direction = 0;

end