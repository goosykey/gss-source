classdef satellite
    %SATELLITE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        orbitepoch_gd;
        starttime_gd;
        endtime_gd;
        ini0; % a e om Om in u
        prop1;
    end
    
    properties (GetAccess = private)
        pro1c
    end
    
    properties (Constant)
        mu = 398600.44;
        R0 = 6378.14;
    end
    
    properties (Dependent)
        ini_cart;
        test_dep;
    end
    
    methods
        
        %% Constructor method
        function obj = satellite (feck) % constructor method
            if nargin > 0
                if isnumeric(feck)
                    obj.prop1 = feck;
                else
                    error('Value must be numeric')
                end
            end
        end
        
        %% Test
        function boob (obj,jojo)
            obj.prop1 = jojo;
            fprintf('fufel\n');
        end
    end
    
end

