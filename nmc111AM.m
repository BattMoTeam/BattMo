classdef nmc111AM < ComponentModel
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
        
        % Physical constants
        con = physicalConstants();
        
        % Lithium data structure
        Li
        
        % Electron data structure
        e

        % Physicochemical properties
        eps         % constant ? 
        spCAh       % Specific Capacity,            [Ah kg^-1]
        rho         % Mass Density,                 [kg m^-3] or [g L^-1]
        theta0      % Minimum lithiation, 0% SOC    [-]
        theta100    % Maximum lithiation, 100% SOC  [-]
        sigma       % Solid conductivity,           [S m^-1]
        cp          % Molar Heat Capacity,          [J kg^-1 K^-1]             
        k0          % Reference rate constant,      [m^2.5 mol^-0.5 s^-1]
        Eak         % Reaction activation energy,   [J mol^-1] 
        Asp         % Surface area,                 [m2 m^-3]
    end
    
    methods

        function model = nmc111AM(name)
        %NMC111 Construct an instance of the nmc111 class
        %   model = nmc111(cLi, T) SOC is the state-of-charge of the
        %   battery (0-1) and T is the temperature in Kelvin [K]  

            model = model@ComponentModel(name);
                
            % primary variables
            model.pnames = {'phi', 'Li'};

            % state variables
            names = {'phi', ...    % Potential
                     'T', ...      % temperature
                     'SOC', ...
                     'Li', ...     % Lithium concentration
                     'OCP', ...    % Open-circuit potential        [V]
                     'k', ...      % Reaction rate constant        [m^2.5 mol^-0.5 s^-1]
                     'D', ...      % Diffusion
                    };
            model.names = names;
            model = model.setupVarDims();
            
            names = {'k', 'D', 'OCP'};
            inputnames = {'T', 'Li'};
            updatefn = @(model, state) model.updateQuantities(state);
            for ind = 1 : numel(names)
                name = names{ind};
                model = model.addPropFunction(name, updatefn, inputnames, {'.'});
            end
            
            
            % Define material constants
            model.spCAh    = 155;      % [Ah kg^-1]
            model.rho      = 4650;     % [kg m^-3]
            model.theta0   = 0.99174;  % at 0% SOC [-]
            model.theta100 = 0.49550;  % at 100% SOC [-]
            model.Li.cmax  = 51554;    % [mol m^-3]
            model.Li.D0    = 1e-14;    % [m^2 s^-1]
            model.Li.EaD   = 5000;     % [J mol^-1]
            model.sigma    = 100;      % [S m^-1]
            model.cp       = 700;      % [J kg^-1 K^-1]
            model.k0       = 2.334e-11;% [m^2.5 mol^-0.5 s^-1]
            model.Eak      = 5000;     % [J mol^-1]
            model.Asp      = 885000;   % [m2 m^-3]
            model.eps      = 0.8;  
            
        end
        
        function state = updateQuantities(model, state)
        %EQUILIBRIUM Calculate the equilibrium properties of the
        %electrode active material
        %   Calculate the equilibrium open cirucuit potential of
        %   nmc111 according to the model used by Torchio et al [1].

            T     = state.T;
            cs    = state.Li;
            
            % Define ideal gas constant
            R = 8.314;      % [J mol^-1 K^-1]
            
            % Define reference temperature
            refT = 298.15;  % [K]
            
            % Calculate reaction rate constant
            k = model.k0 .* exp( -model.Eak ./ model.con.R .* (1./T - 1/refT));
            
            % Calculate solid diffusion coefficient, [m^2 s^-1]
            D = model.Li.D0.*exp(-model.Li.EaD./R*(1./T - 1/refT));
            
            
            % Calculate the lithiation of the active material. This
            % is a simplification for the initial code! The "real"
            % value of theta should be calculated using the surface
            % concentration of Li and the maximum lithium
            % concentration:
            %
            theta = cs ./ model.Li.cmax;
            %model.theta = (model.theta0 - model.theta100) .* model.SOC + model.theta100;
            
            % Calculate the open-circuit potential at the reference
            % temperature for the given lithiation
            refOCP =  ( -4.656 ...
                        + 88.669 .* theta.^2 ...
                        - 401.119 .* theta.^4 ...
                        + 342.909 .* theta.^6 ...
                        - 462.471 .* theta.^8 ...
                        + 433.434 .* theta.^10 ) ./ ...
                      ( -1 ...
                        + 18.933 .* theta.^2 ...
                        - 79.532 .* theta.^4 ...
                        + 37.311 .* theta.^6 ...
                        - 73.083 .* theta.^8 ...
                        + 95.960 .* theta.^10 );

            % Calculate the entropy change at the given lithiation
            dUdT = -1e-3 .* ...
                   (  0.199521039 ...
                      - 0.928373822 .* theta ...
                      + 1.364550689000003 .* theta.^2 ...
                      - 0.611544893999998 .* theta.^3 ) ./ ...
                   (  1 ...
                      - 5.661479886999997 .* theta ...
                      + 11.47636191 .* theta.^2 ...
                      - 9.82431213599998 .* theta.^3 ...
                      + 3.048755063 .* theta.^4 );

            % Calculate the open-circuit potential of the active
            % material
            OCP = refOCP + (T - refT) .* dUdT;

            state.D   = D;
            state.OCP = OCP;
            state.k   = k;
            
        end
        
    end
end

%% References
%   [1] Torchio et al, Journal of The Electrochemical Society, 163 (7)
%   A1192-A1205 (2016), DOI: 10.1149/2.0291607jes

