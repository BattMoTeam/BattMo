classdef PhysicalConstants
    %PHYSICALCONSTANTS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        R       % Ideal gas constant
        F       % Faraday constant
        c0      % Standard concentration
    end
    
    methods
        function obj = PhysicalConstants()
            %PHYSICALCONSTANTS Construct an instance of this class
            %   Detailed explanation goes here
            obj.R = 8.31446261815324;
            obj.F = 96485.3329;
            obj.c0 = 1000;
        end
    end
end

