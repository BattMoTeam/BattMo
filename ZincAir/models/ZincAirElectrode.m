classdef ZincAirElectrode < ElectronicComponent

    properties
        V                    % Molar Volume
        externalCouplingTerm %
    end

    methods

        function model = ZincAirElectrode(inputparams)

            model = model@ElectronicComponent(inputparams);

            fdnames = {'V', ...
                       'externalCouplingTerm'};
            model = dispatchParams(model, inputparams, fdnames);

        end

        function model = registerVarAndPropfuncNames(model)

            %% Declaration of the Dynamical Variables and Function of the model
            % (setup of varnameList and propertyFunctionList)

            model = registerVarAndPropfuncNames@ElectronicComponent(model);

            varnames = {};
            % Volume fraction
            varnames{end + 1} = 'volumeFraction';
            % Accumulation term for active component
            varnames{end + 1} = 'accumTerm';
            % Source term for mass conservation equation
            varnames{end + 1} = 'sourceTerm';
            % Mass conservation equation
            varnames{end + 1} = 'massCons';

            model = model.registerVarNames(varnames);

            fn = @() ZincAirElectrode.updateMassCons;
            inputnames = {'accumTerm', 'sourceTerm'};
            model = model.registerPropFunction({'massCons', fn, inputnames});

            fn = @() ZincAirElectrode.updateConductivity;
            inputnames = {'volumeFraction'};
            model = model.registerPropFunction({'conductivity', fn, inputnames});

            fn = @() ZincAirElectrode.updateAccumTerm;
            fn = {fn, @(propfunc) PropFunction.accumFuncCallSetupFn(propfunc)};
            inputnames = {'volumeFraction'};
            model = model.registerPropFunction({'accumTerm', fn, inputnames});

        end

        function state = updateAccumTerm(model, state, state0, dt)

            vols = model.G.getVolumes();
            state.accumTerm = 1/dt*vols.*(state.volumeFraction - state0.volumeFraction);

        end



        function state = updateConductivity(model, state)

            sigma = model.EffectiveElectricalConductivity; % this prefix "Effective" should be changed in this property name

            vf = state.volumeFraction;

            state.conductivity = (vf.^1.5).*sigma;

        end


        function state = updateMassCons(model, state)
            accumTerm = state.accumTerm;
            sourceTerm = state.sourceTerm;

            state.massCons = accumTerm - sourceTerm;
        end

    end

end



%{
Copyright 2021-2024 SINTEF Industry, Sustainable Energy Technology
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
