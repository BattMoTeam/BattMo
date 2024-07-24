classdef HomogeneousBlock < BaseModel
%
    properties

        ElectronicModel
        ThermalModel        
        Control
        
        couplingTerms % Coupling terms
        use_thermal

        initT % Initial temperature

    end

    methods

        function model = HomogeneousBlock(inputparams)

            model = model@BaseModel();
            fdnames = {'G'                         , ...
                       'couplingTerms'             , ...
                       'initT'                     , ...
                       'use_thermal'};
            
            model.ElectronicModel = ElectronicComponent(inputparams.ElectronicModel);

            if model.use_thermal
                model.ThermalModel = ThermalComponent(inputparams.ThermalModel);
            end

            model.Control = model.setupControl(inputparams.Control);

            model = model.equipModelForComputation();
            
        end

        function model = registerVarAndPropfuncNames(model)

            %% Declaration of the Dynamical Variables and Function of the model
            % (setup of varnameList and propertyFunctionList)

            model = registerVarAndPropfuncNames@BaseModel(model);

            % defines shorthands for the submodels
            el      = 'ElectronicModel';
            thermal = 'ThermalModel';
            ctrl    = 'Control';

            varnames = {{thermal, 'jHeatOhmSource'     }, ...
                        {thermal, 'jHeatChemicalSource'}, ...
                        {thermal, 'jHeatReactionSource'}}

            model = model.removeVarNames(varnames);
            
            fn = @HomogeneousBlock.setupEIEquation;
            inputnames = {{ctrl, 'E'}, ...
                          {ctrl, 'I'}};
            inputnames{end + 1} = {el, 'phi'};
            model = model.registerPropFunction({{ctrl, 'EIequation'}, fn, inputnames});

            fn = @HomogeneousBlock.setupExternalCoupling;
            inputnames = {{el, 'phi'}, ...
                          {el, 'conductivity'}};
            model = model.registerPropFunction({{el, 'jBcSource'}, fn, inputnames});
            model = model.registerPropFunction({{el, 'eSource'}, fn, inputnames}); % set to zero
            if model.use_thermal
                model = model.registerPropFunction({'jFaceBc', fn, inputnames});
            end
            
            
            %% Function that update the Thermal Ohmic Terms

            if model.use_thermal
                
                fn = @HomogeneousBlock.updateJfaceBc;
                varnames ={{el, 'j'}       , ...
                           {el, 'conductivity'}};
                model = model.registerPropFunction({{el, 'jFace'}, fn, inputnames});



                fn = @HomogeneousBlock.updateHeatSourceTerms;
                varnames ={{el, 'jFace'}       , ...
                           {el, 'conductivity'}};
                model = model.registerPropFunction({{thermal, 'jHeatSource'}, fn, inputnames});

            end

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
