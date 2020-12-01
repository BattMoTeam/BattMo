classdef currentCollectorAM < handle

    properties
        % Identification properties
        name        % Name of current collector
        
        % Physical constants
        con = physicalConstants();
        
        % Electron data structure
        e
        
        % State properties
        T           % Temperature,                  [K]
        phi         % Local electric potential,     [V]
        refOCP      % Reference open circuit 
                    % potential at standard
                    % temperature                   [V]
        OCP         % Open-circuit potential        [V]
        
        % State properties (not relevant or not used for the moment)
        dUdT        % Entropy change                [V K^-1]
        SOC         % State of charge               [-]
        theta       % Lithiation                    [-]
        k           % Reaction rate constant        [m^2.5 mol^-0.5 s^-1]
        eps         % Volume fraction,              [-]
        Asp         % Surface area,                 [m2 m^-3]
        
        % Physicochemical properties
        sigma       % Solid conductivity            [S m^-1]
                    
        % Physicochemical properties (not relevant or not used for the moment)
        spCAh       % Specific Capacity             [Ah kg^-1]
        theta0      % Minimum lithiation, 0% SOC    [-]
        theta100    % Maximum lithiation, 100% SOC  [-]
        maxc        % Maximum lithium concentration [mol m^-3]
        rho         % Mass Density                  [kg m^-3] or [g L^-1]
        lambda      % Thermal Conductivity          [W m^-1 K^-1]
        cp          % Molar Heat Capacity           [J kg^-1 K^-1]             
        D0          % Diffusion coefficient         [m^2 s^-1]
        D           % Diffusion coefficient         [m^2 s^-1]
        EaD         % Diffusion activ. energy       [J mol^-1]
        k0          % Reference rate constant       [m^2.5 mol^-0.5 s^-1]
        Eak         % Reaction activation energy    [J mol^-1]
        rp          % Particle radius               [m]
    end
    
    methods
        function obj = currentCollectorAM(T)
            % Construct an instance of the currentCollector class
            %   obj = graphite(T) where T is the temperature in Kelvin [K]
            
            % Set object name
            obj.name = 'current collector';
            
            % Define material constants
            obj.sigma = 500; % [S m^-1]
            obj.T = T;
            
        end
        
    end
end


