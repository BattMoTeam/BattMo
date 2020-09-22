classdef nmc111AM < handle
    %NMC111AM An electrode active material class for electrochemical
    %modelling. 
    %   The nmc111 class describes the properties and
    %   parameterization for the active material of nmc111 electrodes.
    %
    %   The class calculates properties based on experimental
    %   parameterization studies described in the scientific literature.
    %   The validity of the parameterization is limited to the conditions
    %   in which it was reported.
    %
    %   Author: Simon Clark (simon.clark@sintef.no)
    %   Usage:  This code is free to use "as-is" for the purpose of 
    %           research at SINTEF without warranty of any kind. The code
    %           is provided with the hope that it will be helpful. The 
    %           author assumes no liability.         
    %
    %   Acknowledgement: This code builds on the work of many other
    %   scientists over decades of research. Their work is gratefully
    %   acknowledged and cited throughout the code. 
    %
    %   Revision History:
    %       03.06.2020: SC (simon.clark@sintef.no) - New Energy Solutions 
    %                   Initial version (0.0-alpha)
    
    properties
        % Identification properties
        name        % Name of the electrode material
        
        % Physical constants
        con = physicalConstants();
        
        % Lithium data structure
        Li
        
        % Electron data structure
        e
        
        % State properties
        T           % Temperature,                  [K]
        phi         % Local electric potential,     [V]
        refOCP      % Reference open circuit, 
                    % potential at standard,
                    % teperature,                   [V]
        OCP         % Open-circuit potential,       [V]
        dUdT        % Entropy change,               [V K^-1]
        SOC         % State of charge,              [-]
        theta       % Lithiation,                   [-]
        k           % Reaction rate constant,       [m^2.5 mol^-0.5 s^-1]
        eps         % Volume fraction,              [-]
        Asp         % Surface area,                 [m2 m^-3]
        
        % Physicochemical properties
        MW          % Molecular weight,             [kg mol^-1]
        spCAh       % Specific Capacity,            [Ah kg^-1]
        theta0      % Minimum lithiation, 0% SOC    [-]
        theta100    % Maximum lithiation, 100% SOC  [-]
        rho         % Mass Density,                 [kg m^-3] or [g L^-1]
        sigma       % Solid conductivity,           [S m^-1]
        lambda      % Thermal Conductivity,         [W m^-1 K^-1]
        cp          % Molar Heat Capacity,          [J kg^-1 K^-1]             
        k0          % Reference rate constant,      [m^2.5 mol^-0.5 s^-1]
        Eak         % Reaction activation energy,   [J mol^-1]
    end
    
    methods
        function obj = nmc111AM(SOC, T)
            %NMC111 Construct an instance of the nmc111 class
            %   obj = nmc111(cLi, T) SOC is the state-of-charge of the
            %   battery (0-1) and T is the temperature in Kelvin [K]  
            
            % Set object name
            obj.name = 'nmc111';
            
            % Define material constants
            obj.spCAh       = 155;      % [Ah kg^-1]
            obj.rho         = 4650;     % [kg m^-3]
            obj.theta0    = 0.99174;  % at 0% SOC [-]
            obj.theta100    = 0.49550;  % at 100% SOC [-]
            obj.Li.cmax        = 51554;    % [mol m^-3]
            obj.Li.D0          = 1e-14;    % [m^2 s^-1]
            obj.Li.EaD         = 5000;     % [J mol^-1]
            obj.sigma       = 100;      % [S m^-1]
            obj.cp          = 700;      % [J kg^-1 K^-1]
            obj.k0          = 2.334e-11;% [m^2.5 mol^-0.5 s^-1]
            obj.Eak         = 5000;     % [J mol^-1]
            obj.Asp         = 885000;   % [m2 m^-3]
            
            % Define material state
            obj.SOC = SOC;
            m = (1 ./ (obj.theta100 - obj.theta0));
            b = -m .* obj.theta0;
            obj.theta = (obj.SOC - b) ./ m;
            obj.Li.cs = obj.theta .* obj.Li.cmax;
            obj.T = T;      
            
            % Update the properties
            obj.update()
            
            obj.phi = obj.OCP;

        end
        
        function update(obj)
            %UPDATE Update the electrode properties
            %   Calculate the updated properties of the active material
            %   at the given state.
            obj.equilibrium()
            obj.kinetics()
            obj.diffusionCoefficient()
        end
        
        function kinetics(obj)
            %KINETICS Calculate the kinetic parameters for the Li+
            %intercalation reaction
            
            % Define ideal gas constant
            R = 8.314;      % [J mol^-1 K^-1]
            
            % Define reference temperature
            refT = 298.15;  % [K]
            
            % Calculate reaction rate constant
            obj.k = obj.k0 .* exp( -obj.Eak ./ R .* (1./obj.T-1/refT));
            
        end
        
        function diffusionCoefficient(obj)
            %DIFFUSIONCOEFFICIENT Calculate the solid diffusion coefficient of Li+ in
            %the active material
            %   Calculate the solid phase diffusion coefficient of Li+ in
            %   nmc111 according to the model used by Torchio et al [1].
            
            % Define ideal gas constant
            R = 8.314;      % [J mol^-1 K^-1]
            
            % Define reference temperature
            refT = 298.15;  % [K]
           
            % Calculate solid diffusion coefficient, [m^2 s^-1]
            obj.Li.D = obj.Li.D0 .* exp(-obj.Li.EaD./R*(1./obj.T-1/refT));
            
            
        end
        
        function equilibrium(obj)
            %EQUILIBRIUM Calculate the equilibrium properties of the
            %electrode active material
            %   Calculate the equilibrium open cirucuit potential of
            %   nmc111 according to the model used by Torchio et al [1].

            % Set the reference temperature
            refT = 298.15;

                    % Calculate the lithiation of the active material. This
                    % is a simplification for the initial code! The "real"
                    % value of theta should be calculated using the surface
                    % concentration of Li and the maximum lithium
                    % concentration:
                    %
                    obj.theta = obj.Li.cs ./ obj.Li.cmax;
                    %obj.theta = (obj.theta0 - obj.theta100) .* obj.SOC + obj.theta100;
                    
                    % Calculate the open-circuit potential at the reference
                    % temperature for the given lithiation
                    obj.refOCP =  ( -4.656 ...
                        + 88.669 .* obj.theta.^2 ...
                        - 401.119 .* obj.theta.^4 ...
                        + 342.909 .* obj.theta.^6 ...
                        - 462.471 .* obj.theta.^8 ...
                        + 433.434 .* obj.theta.^10 ) ./ ...
                       ( -1 ...
                        + 18.933 .* obj.theta.^2 ...
                        - 79.532 .* obj.theta.^4 ...
                        + 37.311 .* obj.theta.^6 ...
                        - 73.083 .* obj.theta.^8 ...
                        + 95.960 .* obj.theta.^10 );

                    % Calculate the entropy change at the given lithiation
                    obj.dUdT = -1e-3 .* ...
                       (  0.199521039 ...
                        - 0.928373822 .* obj.theta ...
                        + 1.364550689000003 .* obj.theta.^2 ...
                        - 0.611544893999998 .* obj.theta.^3 ) ./ ...
                       (  1 ...
                        - 5.661479886999997 .* obj.theta ...
                        + 11.47636191 .* obj.theta.^2 ...
                        - 9.82431213599998 .* obj.theta.^3 ...
                        + 3.048755063 .* obj.theta.^4 );

                    % Calculate the open-circuit potential of the active
                    % material
                    obj.OCP = obj.refOCP + (obj.T - refT) .* obj.dUdT;

        end
        
    end
end

%% References
%   [1] Torchio et al, Journal of The Electrochemical Society, 163 (7)
%   A1192-A1205 (2016), DOI: 10.1149/2.0291607jes

