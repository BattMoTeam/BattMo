classdef ImpedanceElectrolyser < Electrolyser

    properties

    end

    methods

        function model = ImpedanceElectrolyser(inputparams)

            model = model@Electrolyser(inputparams);

        end


        function model = registerVarAndPropfuncNames(model)

        %% Declaration of the Dynamical Variables and Function of the model
        % (setup of varnameList and propertyFunctionList)

            model = registerVarAndPropfuncNames@Electrolyser(model);

            model = model.registerVarName('omega');

            model = model.setAsStaticVarName('omega');
            
            % defines shorthands for the submodels
            inm = 'IonomerMembrane';
            her = 'HydrogenEvolutionElectrode';
            oer = 'OxygenEvolutionElectrode';

            ctl = 'CatalystLayer';
            exr = 'ExchangeReaction';
            ptl = 'PorousTransportLayer';

            eldes = {her, oer};
            layers = {ctl, exr};
            
            fn = @() ImpedanceElectrolyser.updateIonomerImpedanceTerm;
            inputnames = {{inm, 'H2Oceps'}, 'omega'};
            model = model.registerPropFunction({{inm, 'H2Oaccum'}, fn, inputnames});

            fn = @() ImpedanceElectrolyser.updatePorousTransportLayerimpedanceTerm;
            % HER
            ngas = model.(her).(ptl).gasInd.ngas;
            inputnames = {VarName({her, ptl}, 'compGasMasses', ngas), {her, ptl, 'OHceps'}, 'omega'};
            model = model.registerPropFunction({VarName({her, ptl}, 'compGasAccums', ngas), fn, inputnames});
            model = model.registerPropFunction({{her, ptl, 'OHaccum'}, fn, inputnames});
            % OER
            ngas = model.(her).(ptl).gasInd.ngas;
            inputnames = {VarName({oer, ptl}, 'compGasMasses', ngas), {oer, ptl, 'OHceps'}, 'omega'};
            model = model.registerPropFunction({VarName({oer, ptl}, 'compGasAccums', ngas), fn, inputnames});
            model = model.registerPropFunction({{oer, ptl, 'OHaccum'}, fn, inputnames});
            
            fn = @() ImpedanceElectrolyser.updatePorousTransportLayerLiquidImpedanceTerm;
            % HER
            inputnames = {{her, ptl, 'liqrhoeps'}, 'omega'};
            model = model.registerPropFunction({{her, ptl, 'liquidAccumTerm'}, fn, inputnames});
            % OER
            inputnames = {{oer, ptl, 'liqrhoeps'}, 'omega'};
            model = model.registerPropFunction({{oer, ptl, 'liquidAccumTerm'}, fn, inputnames});

        end


        function state = updateLiquideImpedanceTerm(model, state)

            % NOT UPDATED
            
            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';
            co      = 'Coating';
            am      = 'ActiveMaterial';
            sd      = 'SolidDiffusion';

            op = model.(ne).(co).(am).(sd).operators;

            c     = state.(ne).(co).(am).(sd).c;
            omega = state.Control.omega;
            
            state.(ne).(co).(am).(sd).massAccum = i*omega.*op.vols.*c;

            
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
