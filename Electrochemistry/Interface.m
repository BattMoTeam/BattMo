classdef Interface < BaseModel

    properties

        % Physical constants
        constants = PhysicalConstants();

        % Appelation name of the active material
        name

        cmax
        
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

        updateOCPFunc % Function handler to update OCP
        
    end

    methods

        function model = Interface(paramobj)

            model = model@BaseModel();

             % OBS : All the submodels should have same backend (this is not assigned automaticallly for the moment)
            model.AutoDiffBackend = SparseAutoDiffBackend('useBlocks', false);

            fdnames = {'G'                      , ...
                       'name'                   , ...
                       'specificCapacity'       , ...
                       'rho'                    , ...
                       'theta0'                 , ...
                       'theta100'               , ...
                       'cmax'                     , ...
                       'k0'                     , ...
                       'Eak'                    , ...
                       'rp'                     , ...
                       'volumetricSurfaceArea'  , ...
                       'density'                , ...
                       'n'                      , ...
                       'volumeFraction'};

            model = dispatchParams(model, paramobj, fdnames);

            model.updateOCPFunc = str2func(paramobj.updateOCPFunc.functionname);

        end

        function model = registerVarAndPropfuncNames(model)
            
            %% Declaration of the Dynamical Variables and Function of the model
            % (setup of varnameList and propertyFunctionList)

            model = registerVarAndPropfuncNames@BaseModel(model);
            
            varnames = {};
            % Temperature
            varnames{end + 1} = 'T';
            % Status of Charge
            varnames{end + 1} = 'SOC';
            % potential in electrode
            varnames{end + 1} = 'phiElectrode';
            % charge carrier concentration in electrode - value at surface
            varnames{end + 1} = 'cElectrodeSurface';
            % potential in electrolyte
            varnames{end + 1} = 'phiElectrolyte';
            % charge carrier concentration in electrolyte
            varnames{end + 1} = 'cElectrolyte';
            % eta
            varnames{end + 1} = 'eta';
            % Reaction rate
            varnames{end + 1} = 'R';
            % External potential drop used in Butler-Volmer
            varnames{end + 1} = 'externalPotentialDrop';
            % 
            varnames{end + 1} = 'dUdT';
            % OCP
            varnames{end + 1} = 'OCP';
            % Reaction rate coefficient
            varnames{end + 1} = 'j0';
            
            model = model.registerVarNames(varnames);
            
            fn = @Interface.updateReactionRateCoefficient;
            inputnames = {'T', 'cElectrolyte', 'cElectrodeSurface'};
            model = model.registerPropFunction({'j0', fn, inputnames});

            fn = @Interface.updateOCP;
            inputnames = {'cElectrodeSurface', 'T'};
            model = model.registerPropFunction({'OCP', fn, inputnames});
            % model = model.registerPropFunction({'dUdT', fn, inputnames});
            
            % fn = @Interface.updateEta;
            % inputnames = {'phiElectrolyte', 'phiElectrode', , 'OCP'};
            % model = model.registerPropFunction({'eta', fn, inputnames});
            
            fn = @Interface.updateEtaWithEx;
            inputnames = {'phiElectrolyte', 'phiElectrode', 'OCP', 'externalPotentialDrop'};
            model = model.registerPropFunction({'eta', fn, inputnames});            
            
            fn = @Interface.updateReactionRate;
            inputnames = {'T', 'eta', 'j0'};
            model = model.registerPropFunction({'R', fn, inputnames});
            
            
        end
        
        function state = dipatchTemperature(model, state)

            sd = 'SolidDiffusion';
            state.(sd).T = state.T;
            
        end
        
        function state = updateOCP(model, state)

            c = state.cElectrodeSurface;
            T = state.T;

            cmax = model.cmax;

            func = model.updateOCPFunc;

            [state.OCP, state.dUdT] = func(c, T, cmax);
            
        end

        function state = updateReactionRateCoefficient(model, state)

            Tref = 298.15;  % [K]

            cmax = model.cmax;
            k0   = model.k0;
            Eak  = model.Eak;
            n    = model.n;
            R    = model.constants.R;
            F    = model.constants.F;

            T      = state.T;
            cElyte = state.cElectrolyte;
            c      = state.cElectrodeSurface;
            
            % Calculate reaction rate constant
            k = k0.*exp(-Eak./R.*(1./T - 1/Tref));

            % We use regularizedSqrt to regularize the square root function and avoid the blow-up of derivative at zero.
            th = 1e-3*cmax;
            j0 = k.*regularizedSqrt(cElyte.*(cmax - c).*c, th)*n*F;

            state.j0 = j0;

        end

        function state = updateEta(model, state)

            phiElyte = state.phiElectrolyte;
            phiElde  = state.phiElectrode;
            OCP      = state.OCP;

            state.eta = (phiElde - phiElyte - OCP);

        end
        
        function state = updateEtaWithEx(model, state)

            phiElyte = state.phiElectrolyte;
            phiElde  = state.phiElectrode;
            OCP      = state.OCP;
            dphi     = state.externalPotentialDrop;

            %% FIXME : check sign
            state.eta = (phiElde - phiElyte - OCP - dphi);

        end
        
            
        function state = updateReactionRate(model, state)

            n = model.n;
            F = model.constants.F;

            T   = state.T;
            j0  = state.j0;
            eta = state.eta;
            
            R = model.volumetricSurfaceArea.*ButlerVolmerEquation(j0, 0.5, n, eta, T);

            state.R = R/(n*F); % reaction rate in mole/meter^3/second

        end
        
    end
end

%% References
%   [1] Torchio et al, Journal of The Electrochemical Society, 163 (7)
%   A1192-A1205 (2016), DOI: 10.1149/2.0291607jes



%{
Copyright 2021-2022 SINTEF Industry, Sustainable Energy Technology
and SINTEF Digital, Mathematics & Cybernetics.

This file is part of The Battery Modeling Toolbox BattMo

BattMo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BattMo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BattMo.  If not, see <http://www.gnu.org/licenses/>.
%}
