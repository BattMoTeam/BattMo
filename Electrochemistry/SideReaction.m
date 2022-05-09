classdef SideReaction < BaseModel

    properties

        % Physical constants
        constants = PhysicalConstants();

        beta = 0.5               % side reaction buttler-volmer  coefficient [-]
        k = 1.36e-12             % side reaction rate constant [m/s]
        conductivity = 5e-6;     % ionic conductivity [S/m]
        
    end

    methods

        function model = SideReaction(paramobj)

            model = model@BaseModel();

             % OBS : All the submodels should have same backend (this is not assigned automaticallly for the moment)
            model.AutoDiffBackend = SparseAutoDiffBackend('useBlocks', false);

            fdnames = {'beta', ...
                       'k'   , ...
                       'conductivity'};

            model = dispatchParams(model, paramobj, fdnames);

        end

        function model = registerVarAndPropfuncNames(model)
            
            %% Declaration of the Dynamical Variables and Function of the model
            % (setup of varnameList and propertyFunctionList)

            model = registerVarAndPropfuncNames@BaseModel(model);
            
            varnames = {};
            % potential in electrode
            varnames{end + 1} = 'phiElectrode';
            % charge carrier concentration in electrode - value at surface
            varnames{end + 1} = 'c';
            % potential in electrolyte
            varnames{end + 1} = 'phiElectrolyte';
            % Reaction rate
            varnames{end + 1} = 'R';
            % Reaction rate coefficient
            varnames{end + 1} = 'j0';
            % External potential drop used in Butler-Volmer
            varnames{end + 1} = 'externalPotentialDrop';
            
            model = model.registerVarNames(varnames);
            
            fn = @Interface.updateReactionRateCoefficient;
            inputnames = {'c'};
            model = model.registerPropFunction({'j0', fn, inputnames});

            fn = @Interface.updateReactionRate;
            inputnames = {'phiElectrolyte', 'phiElectrode', 'j0', 'externalPotentialDrop'};
            model = model.registerPropFunction({'R', fn, inputnames});
            
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

        function state = updateReactionRate(model, state)

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
