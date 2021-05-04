classdef ActiveMaterialInputParams 
    
    properties
        
        G

        name
        
        specificCapacity        % [Ah kg^-1]
        rho                     % [kg m^-3]
        theta0                  % at 0% SOC [-]
        theta100                % at 100% SOC[-]
        Li % struct with fields
           %  - cmax            % [mol m^-3]
           %  - D0              % [m^2 s^-1]
           %  - EaD             % [J mol^-1]
        electricalConductivity  % [S m^-1]
        cp                      % [J kg^-1 K^-1]
        k0                      % [m^2.5 mol^-0.5 s^-1]
        Eak                     % [J mol^-1]
        volumetricSurfaceArea   % [m2 m^-3]
        volumeFraction         
        
    end
    
    methods
        
        function paramobj = ActiveMaterialInputParams();
            paramobj.G = struct();
            paramobj.name = char();
            paramobj.Li = struct();
        end
        
    end
    
    
end
