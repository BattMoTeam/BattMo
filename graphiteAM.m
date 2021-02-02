classdef graphiteAM < ComponentModel
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
        
        % Physical constants
        con = physicalConstants();
        
        % Lithium data structure
        Li
        
        % Electron data structure
        e
        
        % Physicochemical properties
        eps
        Asp         % Surface area,                 [m2 m^-3]
        spCAh       % Specific Capacity             [Ah kg^-1]
        theta0      % Minimum lithiation, 0% SOC    [-]
        theta100    % Maximum lithiation, 100% SOC  [-]
        maxc        % Maximum lithium concentration [mol m^-3]
        rho         % Mass Density                  [kg m^-3] or [g L^-1]
        sigma       % Solid conductivity            [S m^-1]
        lambda      % Thermal Conductivity          [W m^-1 K^-1]
        cp          % Molar Heat Capacity           [J kg^-1 K^-1]             
        D0          % Diffusion coefficient         [m^2 s^-1]
        EaD         % Diffusion activ. energy       [J mol^-1]
        k0          % Reference rate constant       [m^2.5 mol^-0.5 s^-1]
        Eak         % Reaction activation energy    [J mol^-1]
        rp          % Particle radius               [m]
        
    end
    
    methods
        function model = graphiteAM(name)
            
        % GRAPHITE Construct an instance of the graphite class
        % model = graphite(SOC, T) SOC is the state of charge of the
        % electrode (0-1) and T is the temperature in Kelvin [K]

            model = model@ComponentModel(name);
                
            % Define material constants
            model.spCAh    = 360;      % [Ah kg^-1]
            model.rho      = 2240;     % [kg m^-3]
            model.theta0   = 0.1429;   % at 0% SOC [-]
            model.theta100 = 0.85510;  % at 100% SOC[-]
            model.Li.cmax  = 30555;    % [mol m^-3]
            model.Li.D0    = 3.9e-14;  % [m^2 s^-1]
            model.Li.EaD   = 5000;     % [J mol^-1]
            model.sigma    = 100;      % [S m^-1]
            model.cp       = 700;      % [J kg^-1 K^-1]
            model.k0       = 5.031e-11;% [m^2.5 mol^-0.5 s^-1]
            model.Eak      = 5000;     % [J mol^-1]
            model.Asp      = 723600;   % [m2 m^-3]
            model.eps      = 0.8;
            
            % primary variables
            names = {'phi', 'Li'};
            model.pnames = names;

            % state variables
            names = {'phi', ...    % Potential
                     'T', ...      % Temperature
                     'SOC', ...
                     'Li', ...     % Lithium concentration
                     'refOCP', ... % Reference open circuit potential at standard  teperature [V]
                     'OCP', ...    % Open-circuit potential        [V]
                     'k', ...      % Reaction rate constant        [m^2.5 mol^-0.5 s^-1]
                     'D', ...      % Diffusion
                     'eps' ...     % Volume fraction,              [-]    
                    };
            model.names = names; 
            
            propfunctions = {};
            names = {'k', 'D', 'refOCP', 'OCP'};
            updatefn = @(model, state) model.updateQuantities(state);
            for ind = 1 : numel(names)
                name = names{ind};
                propfunction = PropFunction(name, updatefn, '.');
                propfunctions{end + 1} = propfunction;
            end

            model.propfunctions = propfunctions;
            
        end
        
        function state = initializeState(model, state)
        % Used only in debugging for the moment

            T = model.getUpdatedProp(state, 'T');
            SOC = model.getUpdatedProp(state, 'SOC');
            
            m     = (1 ./ (model.theta100 - model.theta0));
            b     = -m .* model.theta0;
            theta = (SOC - b) ./ m;
            cs    = theta .* model.Li.cmax;

            state = model.setProp(state, 'Li', cs);
            
            OCP = model.getUpdatedProp(state, 'OCP');
            state = model.setProp(state, 'phi', OCP);
        end

        
        function state = updateQuantities(model, state)
        % Calculate the solid diffusion coefficient of Li+ in the active material
        % Calculate the solid phase diffusion coefficient of Li+ in
        % graphite according to the model used by Torchio et al [1].
            
            % Define reference temperature
            refT = 298.15;  % [K]
            
            T = model.getUpdatedProp(state, 'T');
            
            % Define reference temperature
            refT = 298.15;  % [K]

            % Calculate reaction rate constant
            k = model.k0 .* exp( -model.Eak ./ model.con.R .* (1./T-1/refT));
                
                
            % Calculate solid diffusion coefficient, [m^2 s^-1]
            D = model.Li.D0 .* exp(-model.Li.EaD./model.con.R*(1./T - 1/refT));

            [cs, state]    = model.getUpdatedProp(state, 'Li');
            
            % Set the reference temperature
            refT = 298.15;

            % Calculate the lithiation of the active material. This
            % is a simplification for the initial code! The "real"
            % value of theta should be calculated using the surface
            % concentration of Li and the maximum lithium
            % concentration:
            %
            theta = cs ./ model.Li.cmax;
            
            % Calculate the open-circuit potential at the reference temperature for the given lithiation
            refOCP = (0.7222 ...
                      + 0.1387 .* theta ...
                      + 0.0290 .* theta.^0.5 ...
                      - 0.0172 ./ theta ... 
                      + 0.0019 ./ theta.^1.5 ...
                      + 0.2808 .* exp(0.9-15.*theta) ... 
                      - 0.7984 .* exp(0.4465 .* theta - 0.4108));

            % Calculate the entropy change at the given lithiation
            dUdT = 1e-3 .* ...
                (  0.005269056 ...
                   + 3.299265709 .* theta ...
                   - 91.79325798 .* theta.^2 ...
                   + 1004.911008 .* theta.^3 ...
                   - 5812.278127 .* theta.^4 ...
                   + 19329.75490 .* theta.^5 ...
                   - 37147.89470 .* theta.^6 ...
                   + 38379.18127 .* theta.^7 ...
                   - 16515.05308 .* theta.^8 ) ./ ...
                (  1 ...
                   - 48.09287227 .* theta ...
                   + 1017.234804 .* theta.^2 ...
                   - 10481.80419 .* theta.^3 ...
                   + 59431.30000 .* theta.^4 ...
                   - 195881.6488 .* theta.^5 ...
                   + 374577.3152 .* theta.^6 ...
                   - 385821.1607 .* theta.^7 ...
                   + 165705.8597 .* theta.^8 );

            % Calculate the open-circuit potential of the active material
            OCP = refOCP + (T - refT) .* dUdT;
            
            state = model.setProp(state, 'D', D);
            state = model.setProp(state, 'refOCP', refOCP);
            state = model.setProp(state, 'OCP'   , OCP);
            state = model.setProp(state, 'k', k);
            
        end
        
    end
end

%% References
%   [1] Torchio et al, Journal of The Electrochemical Society, 163 (7)
%   A1192-A1205 (2016), DOI: 10.1149/2.0291607jes

