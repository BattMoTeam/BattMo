classdef graphiteAM < handle
    %GRAPHITEAM An electrode active material class for electrochemical
    %modelling. 
    %   The graphite class describes the properties and
    %   parameterization for the active material of graphite electrodes.
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
    %                   Initial version (0.0.0-alpha)
    
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
        refOCP      % Reference open circuit 
                    % potential at standard
                    % teperature                    [V]
        OCP         % Open-circuit potential        [V]
        dUdT        % Entropy change                [V K^-1]
        SOC         % State of charge               [-]
        theta       % Lithiation                    [-]
        k           % Reaction rate constant        [m^2.5 mol^-0.5 s^-1]
        eps         % Volume fraction,              [-]
        Asp         % Surface area,                 [m2 m^-3]
        
        % Physicochemical properties
        spCAh       % Specific Capacity             [Ah kg^-1]
        theta0      % Minimum lithiation, 0% SOC    [-]
        theta100    % Maximum lithiation, 100% SOC  [-]
        maxc        % Maximum lithium concentration [mol m^-3]
        rho         % Mass Density                  [kg m^-3] or [g L^-1]
        sigma       % Solid conductivity            [S m^-1]
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
        function obj = graphiteAM(SOC, T)
            %GRAPHITE Construct an instance of the graphite class
            %   obj = graphite(SOC, T) SOC is the state of charge of the
            %   electrode (0-1) and T is the temperature in Kelvin [K]
            
            % Set object name
            obj.name = 'graphite';
            
            % Define material constants
            obj.spCAh       = 360;      % [Ah kg^-1]
            obj.rho         = 2240;     % [kg m^-3]
            obj.theta0      = 0.1429;   % at 0% SOC [-]
            obj.theta100    = 0.85510;  % at 100% SOC[-]
            obj.Li.cmax     = 30555;    % [mol m^-3]
            obj.Li.D0       = 3.9e-14;  % [m^2 s^-1]
            obj.Li.EaD      = 5000;     % [J mol^-1]
            obj.sigma       = 100;      % [S m^-1]
            obj.cp          = 700;      % [J kg^-1 K^-1]
            obj.k0          = 5.031e-11;% [m^2.5 mol^-0.5 s^-1]
            obj.Eak         = 5000;     % [J mol^-1]
            obj.Asp         = 723600;   % [m2 m^-3]
            
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
            obj.diffusion()
            obj.kinetics()
        end
        
        function kinetics(obj)
            %KINETICS Calculate the kinetic parameters for the Li+
            %intercalation reaction
            
            % Define reference temperature
            refT = 298.15;  % [K]
            
            % Calculate reaction rate constant
            obj.k = obj.k0 .* exp( -obj.Eak ./ obj.con.R .* (1./obj.T-1/refT));
            
        end
        
        function diffusion(obj)
            %DIFFUSION Calculate the solid diffusion coefficient of Li+ in
            %the active material
            %   Calculate the solid phase diffusion coefficient of Li+ in
            %   graphite according to the model used by Torchio et al [1].
            
            % Define reference temperature
            refT = 298.15;  % [K]
            
            % Calculate solid diffusion coefficient, [m^2 s^-1]
            obj.Li.D = obj.Li.D0 .* exp(-obj.Li.EaD./obj.con.R*(1./obj.T-1/refT));
            
            
        end
        
        function equilibrium(obj)
            %EQUILIBRIUM Calculate the equilibrium properties of the
            %electrode active material
            %   Calculate the equilibrium open cirucuit potential of
            %   graphite according to the model used by Torchio et al [1].

            % Set the reference temperature
            refT = 298.15;

                    % Calculate the lithiation of the active material. This
                    % is a simplification for the initial code! The "real"
                    % value of theta should be calculated using the surface
                    % concentration of Li and the maximum lithium
                    % concentration:
                    %
                    obj.theta = obj.Li.cs ./ obj.Li.cmax;
                    %obj.theta = (obj.theta100 - obj.theta0) .* obj.SOC;
                    
                    % Calculate the open-circuit potential at the reference
                    % temperature for the given lithiation
                    obj.refOCP = 0.7222 ...
                        + 0.1387 .* obj.theta ...
                        + 0.0290 .* obj.theta.^0.5 ...
                        - 0.0172 ./ obj.theta ...
                        + 0.0019 ./ obj.theta.^1.5 ...
                        + 0.2808 .* exp(0.9-15.*obj.theta) ...
                        - 0.7984 .* exp(0.4465 .* obj.theta - 0.4108);

                    % Calculate the entropy change at the given lithiation
                    obj.dUdT = 1e-3 .* ...
                       (  0.005269056 ...
                        + 3.299265709 .* obj.theta ...
                        - 91.79325798 .* obj.theta.^2 ...
                        + 1004.911008 .* obj.theta.^3 ...
                        - 5812.278127 .* obj.theta.^4 ...
                        + 19329.75490 .* obj.theta.^5 ...
                        - 37147.89470 .* obj.theta.^6 ...
                        + 38379.18127 .* obj.theta.^7 ...
                        - 16515.05308 .* obj.theta.^8 ) ./ ...
                       (  1 ...
                        - 48.09287227 .* obj.theta ...
                        + 1017.234804 .* obj.theta.^2 ...
                        - 10481.80419 .* obj.theta.^3 ...
                        + 59431.30000 .* obj.theta.^4 ...
                        - 195881.6488 .* obj.theta.^5 ...
                        + 374577.3152 .* obj.theta.^6 ...
                        - 385821.1607 .* obj.theta.^7 ...
                        + 165705.8597 .* obj.theta.^8 );

                    % Calculate the open-circuit potential of the active
                    % material
                    obj.OCP = obj.refOCP + (obj.T - refT) .* obj.dUdT;
                    
        end
        
    end
end

%% References
%   [1] Torchio et al, Journal of The Electrochemical Society, 163 (7)
%   A1192-A1205 (2016), DOI: 10.1149/2.0291607jes

