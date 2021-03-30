classdef NMC111InputParams < ActiveMaterialInputParams 
    
    
    methods
        
        function paramobj = NMC111InputParams()
            
            paramobj = paramobj@ActiveMaterialInputParams();
            
            % Define material constants
            paramobj.specificCapacity       = 155;      % [Ah kg^-1]
            paramobj.rho                    = 4650;     % [kg m^-3]
            paramobj.theta0                 = 0.99174;  % at 0% SOC [-]
            paramobj.theta100               = 0.49550;  % at 100% SOC [-]
            paramobj.Li.cmax                = 51554;    % [mol m^-3]
            paramobj.Li.D0                  = 1e-14;    % [m^2 s^-1]
            paramobj.Li.EaD                 = 5000;     % [J mol^-1]
            paramobj.electronicConductivity = 100;      % [S m^-1]
            paramobj.cp                     = 700;      % [J kg^-1 K^-1]
            paramobj.k0                     = 2.334e-11;% [m^2.5 mol^-0.5 s^-1]
            paramobj.Eak                    = 5000;     % [J mol^-1]
            paramobj.volumetricSurfaceArea  = 885000;   % [m2 m^-3]
            paramobj.volumeFraction         = 0.8;  
            
            
            paramobj.specificCapacity       = 360;      % [Ah kg^-1]
            paramobj.rho                    = 2240;     % [kg m^-3]
            paramobj.theta0                 = 0.1429;   % at 0% SOC [-]
            paramobj.theta100               = 0.85510;  % at 100% SOC[-]
            paramobj.Li.cmax                = 30555;    % [mol m^-3]
            paramobj.Li.D0                  = 3.9e-14;  % [m^2 s^-1]
            paramobj.Li.EaD                 = 5000;     % [J mol^-1]
            paramobj.electronicConductivity = 100;      % [S m^-1]
            paramobj.cp                     = 700;      % [J kg^-1 K^-1]
            paramobj.k0                     = 5.031e-11;% [m^2.5 mol^-0.5 s^-1]
            paramobj.Eak                    = 5000;     % [J mol^-1]
            paramobj.volumetricSurfaceArea  = 723600;   % [m2 m^-3]
            paramobj.volumeFraction         = 0.8;
            
        end
        
    end
    
end
