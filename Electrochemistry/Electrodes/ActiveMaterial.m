classdef ActiveMaterial < BaseModel

    properties

        % Physical constants
        constants = PhysicalConstants();

        % Appelation name of the active material
        name

        % Lithium data structure
        Li

        % number of electron transfer
        n

        % Physicochemical properties
        volumeFraction
        volumetricSurfaceArea  % Surface area to volume,       [m2 m^-3]
        density                % [kg m^-3]
        theta0                 % Minimum lithiation, 0% SOC    [-]
        theta100               % Maximum lithiation, 100% SOC  [-]
        k0                     % Reference rate constant       [m^2.5 mol^-0.5 s^-1]
        Eak                    % Reaction activation energy    [J mol^-1]
        rp                     % Particle radius               [m]

        updateOCPFunc % Function handler to update OCP
        
        
    end

    methods

        function model = ActiveMaterial(paramobj)

            model = model@BaseModel();

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
                       'density'                , ...
                       'n'                      , ...
                       'volumeFraction'};

            model = dispatchParams(model, paramobj, fdnames);

            model.updateOCPFunc = str2func(paramobj.updateOCPFunc.functionname);

            %% Declaration of the Dynamical Variables and Function of the model
            % (setup of varnameList and propertyFunctionList)

            varnames = {};
            % Temperature
            varnames{end + 1} = 'T';
            % Status of Charge
            varnames{end + 1} = 'SOC';
            % potential in electrode
            varnames{end + 1} = 'phiElectrode';
            % charge carrier concentration in electrode - value at surface
            varnames{end + 1} = 'cElectrode';
            % charge carrier concentration in electrode - Averaged value
            varnames{end + 1} = 'cElectrodeAveraged';
            % potential in electrolyte
            varnames{end + 1} = 'phiElectrolyte';
            % charge carrier concentration in electrolyte
            varnames{end + 1} = 'cElectrolyte';
            % eta
            varnames{end + 1} = 'eta';
            % Reaction rate
            varnames{end + 1} = 'R';
            % Diffusion coefficient
            % Note : This is the inner particle diffusion coefficient
            varnames{end + 1} = 'D';
            % OCP
            varnames{end + 1} = 'OCP';
            % Reaction rate coefficient
            varnames{end + 1} = 'j0';
            % Solid diffusion equation
            varnames{end + 1} = 'solidDiffusionEq';
            
            model = model.registerVarNames(varnames);
            
            fn = @ActiveMaterial.updateReactionRateCoefficient;
            inputnames = {'cElectrode', 'cElectrolyte', 'T'};
            model = model.registerPropFunction({'j0', fn, inputnames});

            fn = @ActiveMaterial.updateDiffusionCoefficient;
            inputnames = {'T'};
            model = model.registerPropFunction({'D', fn, inputnames});

            fn = @ActiveMaterial.updateOCP;
            inputnames = {'cElectrode', 'T'};
            model = model.registerPropFunction({'OCP', fn, inputnames});

            fn = @ActiveMaterial.updateReactionRate;
            inputnames = {'T', 'phiElectrolyte', 'phiElectrode', 'cElectrode', 'cElectrolyte', 'OCP', 'k'};
            model = model.registerPropFunction({'R', fn, inputnames});
            % model = model.registerPropFunction({'eta', fn, inputnames});
            
            fn = @ActiveMaterial.assembleSolidDiffusionEquation;
            inputnames = {'D', 'R', 'cElectrode', 'cElectrodeAveraged'};
            model = model.registerPropFunction({'solidDiffusionEq', fn, inputnames});

        end

        function state = updateOCP(model, state)
            c = state.cElectrode;
            T = state.T;

            cmax = model.Li.cmax;

            func = model.updateOCPFunc;

            [state.OCP, state.dUdT] = func(c, T, cmax);
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
        % The surface concentration value is computed following polynomial method, as described in ref1 (see below)

            csurf = state.cElectrode;
            cavg = state.cElectrodeAveraged;
            D = state.D;
            R = state.R;

            rp = model.rp;
            a = model.volumetricSurfaceArea;
            % We divide with volumetricSurfaceArea because it was added in the definition of R 
            state.solidDiffusionEq = csurf - cavg + (rp.*R)./(5*a*D);

        end

    end
end

%% References
%   [1] Torchio et al, Journal of The Electrochemical Society, 163 (7)
%   A1192-A1205 (2016), DOI: 10.1149/2.0291607jes




%{
Copyright 2009-2021 SINTEF Industry, Sustainable Energy Technology
and SINTEF Digital, Mathematics & Cybernetics.

This file is part of The Battery Modeling Toolbox BatMo

BatMo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BatMo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BatMo.  If not, see <http://www.gnu.org/licenses/>.
%}
