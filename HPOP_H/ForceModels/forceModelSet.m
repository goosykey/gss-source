classdef forceModelSet < handle
    %FORCEMODELSET Force models for HPOP earth324
    %   Detailed explanation goes here
    
    properties
        grav, maxdeg, maxord, atm, sunmoon, tides, emp
    end
    
    methods
        function obj = forceModelSet(varargin)
            % this constructor method sets default values for force models and
            % calls setter methods for selected inputs
            
            obj.grav = 'full';
            obj.maxdeg = 10;
            obj.maxord = 10;
            obj.atm = 'GOST';
            obj.sunmoon = 'full';
            obj.tides = 'full';
            obj.emp = 'full';
            
            while ~isempty(varargin)
                switch lower(varargin{1})
                    case {'gravity', 'grav', 'nonspherics'}          
                        obj.grav = varargin{2};
                    case {'maxdeg', 'degree'}
                        mustBeNumeric(varargin{2});
                        obj.maxdeg = varargin{2};                        
                    case {'maxord', 'order'}
                        obj.maxord = varargin{2};
                    case {'atmosphere', 'atm', 'drag'}
                        obj.atm = varargin{2};
                    case {'sunmoon', '3body', '3rdbody'}
                        obj.sunmoon = varargin{2};
                    case {'tides', 'tidal'}
                        obj.tides = varargin{2};
                    case {'emp', 'empressure'}
                        obj.emp = varargin{2};
                    otherwise
                        error(['Unexpected option (forceModelSet parser): ' varargin{1}])
                end
                varargin(1:2) = [];
            end
            
        end
        
        
        % BELOW: setter methods with checking argument validity
        
        function set.grav(obj,val)
            mustBeMember(lower(val),{'full','j2','none'});
            obj.grav = val;
        end
        
        function set.maxdeg(obj,val)
            mustBeNumeric(val);
            if mod(val,1) == 0 && val >= 0 && val <= 10
                obj.maxdeg = val;
            else
                error(['Invalid maxdeg value: ' val '; only integers from 0 to 10 are supported'])
            end
        end
        
        function set.maxord(obj,val)
            mustBeNumeric(val);
            if mod(val,1) == 0 && val >= 0 && val <= 10
                obj.maxord = val;
            else
                error(['Invalid maxord value: ' val '; only integers from 0 to 10 are supported'])
            end
        end
        
        function set.atm(obj,val)
            mustBeMember(lower(val),{'gost','j77','jr','none'});
            obj.atm = val;
        end
        
        function set.sunmoon(obj,val)
            mustBeMember(lower(val),{'sun','moon','full','none'});
            obj.sunmoon = val;
        end
        
        function set.tides(obj,val)
            mustBeMember(lower(val),{'full','none'});
            obj.tides = val;
        end
        
        function set.emp(obj,val)
            mustBeMember(lower(val),{'full','none'});
            obj.emp = val;
        end
        
    end
end