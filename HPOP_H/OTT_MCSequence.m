classdef OTT_MCSequence < handle
    %OTTMCS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        SEGMENTS, propTimes, stopConditions, ephBySegment, h = 10, jtBySegment,...
            hasRun = false
    end
    
    methods
        
        %% CONSTRUCT
        
        function obj = OTT_MCSequence()
            obj.propTimes = [];
            obj.stopConditions = {};
            
            obj.ephBySegment = {};
            
            obj.addSegment('type','state');
            obj.addSegment('type','prop');
        end
        
        %% ADD SEGMENT
        
        function addSegment(obj,varargin)
            segment = OTTSegment(varargin{:});
            obj.SEGMENTS = [obj.SEGMENTS,segment];
            currentSegmentNo = numel(obj.SEGMENTS);
            switch lower(segment.type)
                case 'state'
                    obj.propTimes(currentSegmentNo) = nan;
                    obj.stopConditions{currentSegmentNo} = nan;
                case 'prop'
                    obj.propTimes(currentSegmentNo) = 2000;
                    obj.stopConditions{currentSegmentNo} = {'duration', 2000};
                case 'iburn'
                    obj.propTimes(currentSegmentNo) = nan;
                    obj.stopConditions{currentSegmentNo} = nan;
                case 'fburn'
                    obj.propTimes(currentSegmentNo) = 2000;
                    obj.stopConditions{currentSegmentNo} = {'duration', 2000};
                    obj.SEGMENTS(currentSegmentNo).TProp_PF = ...
                        {'f_thrustflag', @(t,y)thrustflag_increase_inclination(t,y), ...
                        'f_thrustvector', @(t,y)thrustVector_increase_inclination(t,y)};
            end
        end
        
        %% RUN NOMINAL SEQUENCE
        
        function run(obj)
            numSeg = numel(obj.SEGMENTS);
            currentState = [];
            currentjd = 2451545;
            for i = 1:numSeg
                seg = obj.SEGMENTS(i);
                
                switch lower(seg.type)
                    
                    case 'state'
                        if seg.TState_ifKepl
                            kep = seg.TState_ini0;
                            [r,v] = math.randv(kep(1),kep(2),kep(5),kep(4),kep(3),kep(6)-kep(3));
                            currentState = [r',v',kep(7)];
                        else
                            currentState = seg.TState_ini0;
                        end
                        currentjd = seg.TState_jd;
                        obj.ephBySegment{i} = currentState;
                        obj.jtBySegment{i} = currentjd;
                        
                    case 'prop'
                        options = odeset('MaxStep', 30,'RelTol', 1e-05);
                        parameters = [seg.TProp_PF,'jd0',currentjd,'forcemodels',{seg.TProp_FM}];
                        tspan = 0:obj.h:1e05;
                        
                        switch obj.stopConditions{i}{1}
                            case 'duration'
                                tspan = 0:obj.h:obj.stopConditions{i}{2};
                            case 'apogee'
                                options.Events = @event_apogee;
                                %add other stop condition cases here
                        end
                        
                        [ti, yi] = ode45( @(t,rv)earth324(t,rv,parameters{:}),...
                            tspan, currentState',options);
                        jti = ti/86400 + currentjd;
                        obj.ephBySegment{i} = yi;
                        obj.jtBySegment{i} = jti;
                        currentjd = jti(end);
                        currentState = yi(end,:);
                        
                    case 'iburn'
                        %TiBurn_ifTSW, TiBurn_dV,
                        r = currentState(1:3);
                        v = currentState(4:6);
                        
                        if seg.TiBurn_ifTSW
                            [~, ~, ~, Om, in, u] = math.cart2kep(r(1), r(2), r(3), v(1), v(2), v(3));
                            v_tsw = math.xyz2tsw(v',Om, in, u) + seg.TiBurn_dV; % ! column vectors here 
                            v = math.tsw2xyz(v_tsw,Om,in,u);
                        else
                            v = v' + seg.TiBurn_dV;
                        end
                        currentState(1:6) = [r, v']; % add mass expulsion effects !
                        obj.ephBySegment{i} = currentState;
                        obj.jtBySegment{i} = currentjd;
                        
                    case 'fburn'
                        options = odeset('MaxStep', 30,'RelTol', 1e-05);
                        parameters = [seg.TProp_PF,'jd0',currentjd,'forcemodels',{seg.TProp_FM}];
                        tspan = 0:obj.h:1e05;
                        
                        switch obj.stopConditions{i}{1}
                            case 'duration'
                                tspan = 0:obj.h:obj.stopConditions{i}{2};
                            case 'apogee'
                                options.Events = @event_apogee;
                                %add other stop condition cases here
                        end
                        
                        [ti, yi] = ode45( @(t,rv)earth324(t,rv,parameters{:}),...
                            tspan, currentState',options);
                        jti = ti/86400 + currentjd;
                        obj.ephBySegment{i} = yi;
                        obj.jtBySegment{i} = jti;
                        currentjd = jti(end);
                        currentState = yi(end,:);
                        
                end
                
            end
            
            obj.hasRun = true;
            
        end
        
        %% PLOT TRAJECTORY
        
        function plotTrajectory(obj)
            vis_earthdraw(obj.jtBySegment{end}(end));
            numSeg = numel(obj.SEGMENTS);
            for i = 1:numSeg
                plot3(...
                    obj.ephBySegment{i}(:,1),...
                    obj.ephBySegment{i}(:,2),...
                    obj.ephBySegment{i}(:,3),...
                    obj.SEGMENTS(i).linecolor, 'linewidth', 2);
            end
            
        end
        
        %% CLEAR EPHEMERIS
        
        function ephClear(obj)
            obj.ephBySegment = {};
        end
        
    end % end methods
end % end classdef
