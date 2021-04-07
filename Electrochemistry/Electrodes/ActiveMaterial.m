classdef ActiveMaterial < PhysicalModel
    
    properties
        
        % Physical constants
        constants = PhysicalConstants();

        % Appelation name of the active material
        name 
        
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
        
    end
    
    methods

        function model = ActiveMaterial(paramobj)
            
            model = model@PhysicalModel([]);
            
             % OBS : All the submodels should have same backend (this is not assigned automaticallly for the moment)
            model.AutoDiffBackend = SparseAutoDiffBackend('useBlocks', false);
            
            fdnames = {'G', ...
                       'name', ...
                       'specificCapacity', ...
                       'rho', ...
                       'theta0', ...
                       'theta100', ...
                       'Li', ...
                       'electronicConductivity', ...
                       'cp', ...
                       'k0', ...
                       'Eak', ...
                       'volumetricSurfaceArea', ...
                       'volumeFraction'};
            
            model = dispatchParams(model, paramobj, fdnames);
            
        end
        
        function state = updateDiffusionConductivityCoefficients(model, state)

            % Define reference temperature
            refT = 298.15;  % [K]

            T = state.T;

            R = model.constants.R;
            
            % Calculate reaction rate constant
            k = model.k0.*exp(-model.Eak./R .*(1./T - 1/refT));
                
            % Calculate solid diffusion coefficient, [m^2 s^-1]
            D = model.Li.D0.*exp(-model.Li.EaD./R*(1./T - 1/refT));

            state.k = k;
            state.D = D;
            
        end
        
        function state = updateReactionRate(model, state);
            
            T        = state.T;
            phiElyte = state.phiElectrolyte;
            % cElyte = state.cElectrolyte; % not used for the moment
            phi      = state.phi;
            % c      = state.c; % not used for the moment
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

