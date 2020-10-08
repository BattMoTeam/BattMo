classdef celgard2500 < FvModel
    %UNTITLED10 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        % Physical constants
        con = physicalConstants()
        
        % Physicochemical properties
        t       % Thickness,        [m]
        eps     % Volume fraction,  [-]
        void    % Porosity,         [-]
        rp      % Pore radius,      [m]
        G       % Gurley number,    [s]
        
        % Mesh properties
        X
        Xb
        N
        dombin
     
        
    end
    
    methods
        function obj = celgard2500()
            %UNTITLED10 Construct an instance of this class
            %   Detailed explanation goes here
            obj.t       = 25e-6;
            obj.void    = 0.55;
            obj.eps     = 1 - obj.void;
            obj.rp      = 0.064e-6 ./ 2;
            obj.G       = 200;
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

