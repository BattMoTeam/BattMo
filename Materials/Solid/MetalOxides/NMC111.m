classdef NMC111 < PhysicalModel
    %NMC111 An electrode active material class for electrochemical
    %modelling. 
    %   The NMC111 class describes the properties and
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
        con = PhysicalConstants();
        
        % Lithium data structure
        Li
        
        % Electron data structure
        e

        % Physicochemical properties
        volumeFraction         % constant ? 
        specificCapacity       % Specific Capacity,            [Ah kg^-1]
        rho         % Mass Density,                 [kg m^-3] or [g L^-1]
        theta0      % Minimum lithiation, 0% SOC    [-]
        theta100    % Maximum lithiation, 100% SOC  [-]
        electronicConductivity       % Solid conductivity,           [S m^-1]
        cp          % Molar Heat Capacity,          [J kg^-1 K^-1]             
        k0          % Reference rate constant,      [m^2.5 mol^-0.5 s^-1]
        Eak         % Reaction activation energy,   [J mol^-1] 

        volumetricSurfaceArea         % Surface area,                 [m2 m^-3]
    end
    
    methods

        function model = NMC111()
        %NMC111 Construct an instance of the nmc111 class
        %   model = nmc111(cLi, T) SOC is the state-of-charge of the
        %   battery (0-1) and T is the temperature in Kelvin [K]  

            model = model@PhysicalModel([]);
            
            % Define material constants
            model.specificCapacity          = 155;      % [Ah kg^-1]
            model.rho                       = 4650;     % [kg m^-3]
            model.theta0                    = 0.99174;  % at 0% SOC [-]
            model.theta100                  = 0.49550;  % at 100% SOC [-]
            model.Li.cmax                   = 51554;    % [mol m^-3]
            model.Li.D0                     = 1e-14;    % [m^2 s^-1]
            model.Li.EaD                    = 5000;     % [J mol^-1]
            model.electronicConductivity    = 100;      % [S m^-1]
            model.cp                        = 700;      % [J kg^-1 K^-1]
            model.k0                        = 2.334e-11;% [m^2.5 mol^-0.5 s^-1]
            model.Eak                       = 5000;     % [J mol^-1]
            model.volumetricSurfaceArea     = 885000;   % [m2 m^-3]
            model.volumeFraction            = 0.8;  
            
        end
        
        function state = updateQuantities(model, state)
        %EQUILIBRIUM Calculate the equilibrium properties of the
        %electrode active material
        %   Calculate the equilibrium open cirucuit potential of
        %   nmc111 according to the model used by Torchio et al [1].

            T  = state.T;
            cs = state.Li;
            
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
            state.theta=theta;
            %model.theta = (model.theta0 - model.theta100) .* model.SOC + model.theta100;
            
            % Calculate the open-circuit potential at the reference
            % temperature for the given lithiation

                    
            coeff1_refOCP = [ -4.656   , ...
                              0        , ...
                              + 88.669 , ...
                              0        , ...
                              - 401.119, ...
                              0        , ...
                              + 342.909, ...
                              0        , ...
                              - 462.471, ...
                              0        , ...
                              + 433.434];
            
            coeff2_refOCP =[ -1      , ...
                             0       , ...
                             + 18.933, ...
                             0       , ...
                             - 79.532, ...
                             0       , ...
                             + 37.311, ...
                             0       , ...
                             - 73.083, ...
                             0       , ...
                             + 95.960];
            
            refOCP = polyval(coeff1_refOCP(end:-1:1),theta)./ polyval(coeff2_refOCP(end:-1:1),theta);    
                    
            % Calculate the entropy change at the given lithiation
            
            coeff1_dUdT = [0.199521039        , ...
                           - 0.928373822      , ...
                           + 1.364550689000003, ...
                           - 0.611544893999998];
            
            coeff2_dUdT = [1                  , ...
                           - 5.661479886999997, ...
                           + 11.47636191      , ... 
                           - 9.82431213599998 , ...
                           + 3.048755063];
            
            dUdT = -1e-3.*polyval(coeff1_dUdT(end:-1:1),theta)./ polyval(coeff2_dUdT(end:-1:1),theta);
            
            % Calculate the open-circuit potential of the active material
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

