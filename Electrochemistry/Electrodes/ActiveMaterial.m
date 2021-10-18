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

        % number of electron transfer
        n

        % Physicochemical properties
        volumeFraction
        volumetricSurfaceArea  % Surface area,                 [m2 m^-3]
        specificCapacity       % Specific Capacity             [Ah kg^-1]
        theta0                 % Minimum lithiation, 0% SOC    [-]
        theta100               % Maximum lithiation, 100% SOC  [-]
        rho                    % Mass Density                  [kg m^-3] or [g L^-1]
        cp                     % Molar Heat Capacity           [J kg^-1 K^-1]
        k0                     % Reference rate constant       [m^2.5 mol^-0.5 s^-1]
        Eak                    % Reaction activation energy    [J mol^-1]
        rp                     % Particle radius               [m]

        updateOCPFunc % Function handler to update OCP
    end

    methods

        function model = ActiveMaterial(paramobj)

            model = model@PhysicalModel([]);

             % OBS : All the submodels should have same backend (this is not assigned automaticallly for the moment)
            model.AutoDiffBackend = SparseAutoDiffBackend('useBlocks', false);

            fdnames = {'G'                      , ...
                       'name'                   , ...
                       'specificCapacity'       , ...
                       'rho'                    , ...
                       'theta0'                 , ...
                       'theta100'               , ...
                       'Li'                     , ...
                       'cp'                     , ...
                       'k0'                     , ...
                       'Eak'                    , ...
                       'rp'                     , ...
                       'volumetricSurfaceArea'  , ...
                       'n'                      , ...
                       'volumeFraction'};

            model = dispatchParams(model, paramobj, fdnames);

            model.updateOCPFunc = str2func(paramobj.updateOCPFunc.functionname);

        end

        function state = updateOCP(model, state)
            c = state.cElectrode;
            T = state.T;

            cmax = model.Li.cmax;

            func = model.updateOCPFunc;

            state.OCP = func(c, T, cmax);
        end

        function state = updateReactionRateCoefficient(model, state)

            Tref = 298.15;  % [K]

            cmax = model.Li.cmax;
            k0   = model.k0;
            Eak  = model.Eak;
            n    = model.n;
            R    = model.constants.R;
            F    = model.constants.F;

            T      = state.T;
            cElyte = state.cElectrolyte;
            c      = state.cElectrode;
            
            % Calculate reaction rate constant
            k = k0.*exp(-Eak./R.*(1./T - 1/Tref));

            % We use regularizedSqrt to regularize the square root function and avoid the blow-up of derivative at zero.
            th = 1e-3*cmax;
            j0 = k.*regularizedSqrt(cElyte.*(cmax - c).*c, th)*n*F;

            state.j0 = j0;

        end

        function state = updateDiffusionCoefficient(model, state)

            Tref = 298.15;  % [K]

            T = state.T;

            R   = model.constants.R;
            D0  = model.Li.D0;
            EaD = model.Li.EaD;

            % Calculate solid diffusion coefficient, [m^2 s^-1]
            D = D0.*exp(-EaD./R*(1./T - 1/Tref));

            state.D = D;
        end

        function state = updateReactionRate(model, state);

            n = model.n;
            F = model.constants.F;

            T        = state.T;
            phiElyte = state.phiElectrolyte;
            phiElde  = state.phiElectrode;
            OCP      = state.OCP;
            j0       = state.j0;

            eta = (phiElde - phiElyte - OCP);
            state.eta = eta;

            R = model.volumetricSurfaceArea.*ButlerVolmerEquation(j0, 0.5, n, eta, T);

            state.R = R/(n*F); % reaction rate in mole/meter^3/second

        end

        function state = assembleSolidDiffusionEquation(model, state)
        % We update the surface concentration of the charge carrier in the active material.
        % The surface concentration value is computed following polynomial method, as described in ref1 (see header)

            csurf = state.cElectrode;
            cavg = state.cElectrodeAveraged;
            D = state.D;
            R = state.R;

            rp = model.rp;
            a = model.volumetricSurfaceArea;

            state.solidDiffusionEq = csurf - cavg + (rp.*R)./(5*a*D);

        end

    end
end

%% References
%   [1] Torchio et al, Journal of The Electrochemical Society, 163 (7)
%   A1192-A1205 (2016), DOI: 10.1149/2.0291607jes

