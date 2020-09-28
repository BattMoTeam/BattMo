classdef seiAM
    %SEIAM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        % Physical constants
        con = physicalConstants();
        
        % State properties
        t       % Thickness,        [m]
        eps     % Volume fraction,  [-]
        
        % Physicochemical properties
                
        
    end
    
    methods
        function obj = seiAM()
            %UNTITLED8 Construct an instance of this class
            %   Detailed explanation goes here
            obj.eps = 0;
            obj.t = 0;
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

