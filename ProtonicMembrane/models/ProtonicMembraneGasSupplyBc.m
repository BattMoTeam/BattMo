classdef ProtonicMembraneGasSupplyBc < BaseModel
    
    properties

        molecularWeights
        T
        
        nGas   % Number of gas (each of them will have a partial pressure). Only needed when gasSupplyType == 'coupled'
        gasInd % Structure whose fieldname give index number of the corresponding gas component.

        constants
        
    end
    
    methods
        
        function model = ProtonicMembraneGasSupplyBc(inputparams)
            
            model = model@BaseModel();
            
            fdnames = {'molecularWeights', ...
                       'T'};

            model = dispatchParams(model, inputparams, fdnames);
            
            model.gasInd.H2O = 1;
            model.gasInd.O2  = 2;
            model.nGas = 2;

            model.constants = PhysicalConstants();
            
        end

        function model = registerVarAndPropfuncNames(model)
            
            model = registerVarAndPropfuncNames@BaseModel(model);

            nGas = model.nGas;
            
            varnames = {};

            % Volume fractions
            varnames{end + 1} = VarName({}, 'massfractions', nGas);
            % Total pressure
            varnames{end + 1} = 'pressure';
            % Density
            varnames{end + 1} = 'density';
            % Densities
            varnames{end + 1} = VarName({}, 'densities', nGas);
            % Mass Flux terms at boundary (outward)
            varnames{end + 1} = VarName({}, 'massFluxes', nGas);
            % Boundary Equation relating boundary massFluxes and boundary pressure values
            varnames{end + 1} = VarName({}, 'bcFluxEquations', nGas);
            % Partial pressures
            varnames{end + 1} = VarName({}, 'pressures', nGas);
            
            model = model.registerVarNames(varnames);

            model = model.setAsExtraVarName(VarName({}, 'pressures', nGas));
            
            fn = @ProtonicMembraneGasSupply.updateMassFraction;
            inputvarnames = {VarName({}, 'massfractions', nGas, 1)};
            outputvarname = VarName({}, 'massfractions', nGas, 2);
            model = model.registerPropFunction({outputvarname, fn, inputvarnames});
                
            fn = @(model, state) ProtonicMembraneGasSupply.updateDensity(model, state);
            fn = {fn, @(prop) PropFunction.literalFunctionCallSetupFn(prop)};
            inputvarnames = {VarName({}, 'massfractions', nGas), 'pressure'};
            model = model.registerPropFunction({'density', fn, inputvarnames});
            
            fn = @(model, state) ProtonicMembraneGasSupply.updateDensities(model, state);
            fn = {fn, @(prop) PropFunction.literalFunctionCallSetupFn(prop)};
            for igas = 1 : nGas
                inputvarnames = {'density', ...
                                 VarName({}, 'massfractions', nGas, igas)};
                outputvarname = VarName({}, 'densities', nGas, igas);
                model = model.registerPropFunction({outputvarname, fn, inputvarnames});
            end

            fn = @(model, state) ProtonicMembraneGasSupply.updatePressures(model, state);
            fn = {fn, @(prop) PropFunction.literalFunctionCallSetupFn(prop)};
            inputvarnames = {VarName({}, 'massfractions', nGas), ...
                             'pressure'};
            outputvarname = VarName({}, 'pressures', nGas);
            model = model.registerPropFunction({outputvarname, fn, inputvarnames});

        end
        
    end

    methods


        function state = updateMassFraction(model, state)

            state.massfractions{2} = 1 - state.massfractions{1};
            
        end

    end
    
end

%{
Copyright 2021-2026 SINTEF Industry, Sustainable Energy Technology
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
