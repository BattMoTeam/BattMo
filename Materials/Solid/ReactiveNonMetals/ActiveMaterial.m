classdef ActiveMaterial < PhysicalModel
    
    properties
        
        % Physical constants
        constants = PhysicalConstants();
        
        % Lithium data structure
        Li
        
        % Electron data structure
        e
        
        % Physicochemical properties
        volumeFraction
        volumetricSurfaceArea  % Surface area,                 [m2 m^-3]
        specificCapacity       % Specific Capacity             [Ah kg^-1]
        theta0                 % Minimum lithiation, 0% SOC    [-]
        theta100               % Maximum lithiation, 100% SOC  [-]
        maxc                   % Maximum lithium concentration [mol m^-3]
        rho                    % Mass Density                  [kg m^-3] or [g L^-1]
        electronicConductivity % Solid conductivity            [S m^-1]
        lambda                 % Thermal Conductivity          [W m^-1 K^-1]
        cp                     % Molar Heat Capacity           [J kg^-1 K^-1]             
        D0                     % Diffusion coefficient         [m^2 s^-1]
        EaD                    % Diffusion activ. energy       [J mol^-1]
        k0                     % Reference rate constant       [m^2.5 mol^-0.5 s^-1]
        Eak                    % Reaction activation energy    [J mol^-1]
        rp                     % Particle radius               [m]
        
        % Names for book-keeping (may not be used)
        ionName
        
    end
    
    methods

        function model = ActiveMaterial(paramobj)
            
            model = model@PhysicalModel([]);
            
            fdnames = {'G'};
            model = dispatchParams(model, paramobj, fdnames);
            
        end
        
        function state = updateMaterialProperties(model, state)
        % Calculate the solid diffusion coefficient of Li+ in the active material
        % Calculate the solid phase diffusion coefficient of Li+ in
        % graphite according to the model used by Torchio et al [1].
            
            % Define reference temperature
            refT = 298.15;  % [K]
            
            T = state.T;
            
            % Define reference temperature
            refT = 298.15;  % [K]

            % Calculate reaction rate constant
            k = model.k0 .* exp( -model.Eak ./ model.con.R .* (1./T-1/refT));
                
                
            % Calculate solid diffusion coefficient, [m^2 s^-1]
            D = model.Li.D0 .* exp(-model.Li.EaD./model.con.R*(1./T - 1/refT));

            cs = state.Li;
            
            % Set the reference temperature
            refT = 298.15;

            % Calculate the lithiation of the active material. This
            % is a simplification for the initial code! The "real"
            % value of theta should be calculated using the surface
            % concentration of Li and the maximum lithium
            % concentration:
            %
            theta = cs ./ model.Li.cmax;
            state.theta = theta;
            % Calculate the open-circuit potential at the reference temperature for the given lithiation
            refOCP = (0.7222 ...
                      + 0.1387 .* theta ...
                      + 0.0290 .* theta.^0.5 ...
                      - 0.0172 ./ theta ... 
                      + 0.0019 ./ theta.^1.5 ...
                      + 0.2808 .* exp(0.9-15.*theta) ... 
                      - 0.7984 .* exp(0.4465 .* theta - 0.4108));
                  
            coeff1 = [0.005269056 ,...
                      + 3.299265709,...
                      - 91.79325798,...
                      + 1004.911008,...
                      - 5812.278127,...
                      + 19329.75490,...
                      - 37147.89470,...
                      + 38379.18127,...
                      - 16515.05308];
               
            coeff2= [1, ...
                     - 48.09287227,...
                     + 1017.234804,...
                     - 10481.80419,...
                     + 59431.30000,...
                     - 195881.6488,...
                     + 374577.3152,...
                     - 385821.1607,...
                     + 165705.8597];
            
            dUdT = 1e-3.*polyval(coeff1(end:-1:1),theta)./ polyval(coeff2(end:-1:1),theta);

            % Calculate the open-circuit potential of the active material
            OCP = refOCP + (T - refT) .* dUdT;
            
            state.D = D;
            state.OCP = OCP;
            state.k = k;
            
        end
        
        function state = updateReactionRate(model, state);
            
            T        = state.T;
            phiElyte = state.phiElectrolyte;
            csElyte  = state.csElectrolyte; % not used for the moment
            phi      = state.phi;
            cs       = state.cs;
            OCP      = state.OCP;
            k        = state.k;
            
            eta = (phi - phiElyte - OCP);
            F = model.constants.F;
            R = model.volumetricSurfaceArea.*ButlerVolmerEquation(k.*model.constants.F, 0.5, 1, eta, T)/F;
             
            state.R = R;
        end
        
        
    end
end

%% References
%   [1] Torchio et al, Journal of The Electrochemical Society, 163 (7)
%   A1192-A1205 (2016), DOI: 10.1149/2.0291607jes

