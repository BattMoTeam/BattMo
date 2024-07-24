classdef HomogeneousBlock < BaseModel
%
    properties

        ElectronicModel
        ThermalModel        
        Control
        
        couplingTerms % Coupling terms
        use_thermal

        initT % Initial temperature

        couplingInd
        
    end

    methods

        function model = HomogeneousBlock(inputparams)

            model = model@BaseModel();
            fdnames = {'couplingTerms'             , ...
                       'initT'                     , ...
                       'use_thermal'};

            model = dispatchParams(model, inputparams, fdnames);
            
            model.ElectronicModel = ElectronicComponent(inputparams.ElectronicModel);

            if model.use_thermal
                model.ThermalModel = ThermalComponent(inputparams.ThermalModel);
            end

            model.Control = HomogeneousBlockControlModel(inputparams.Control);

            coupnames = cellfun(@(coupterm) coupterm.name, model.couplingTerms, 'uniformoutput', false);
            
            couplingInd.input  = find(strcmp(coupnames, 'input'));
            couplingInd.output = find(strcmp(coupnames, 'output'));

            model.couplingInd = couplingInd;
            
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

            if model.use_thermal

                varnames = {{thermal, 'jHeatOhmSource'     }, ...
                            {thermal, 'jHeatChemicalSource'}, ...
                            {thermal, 'jHeatReactionSource'}};

                model = model.removeVarNames(varnames);
                
            else
                
                model = model.removeVarName({el, 'T'});
                
            end

            fn = @HomogeneousBlock.setupEIequation;
            inputnames = {{ctrl, 'E'}, ...
                          {ctrl, 'I'}};
            inputnames{end + 1} = {el, 'phi'};
            model = model.registerPropFunction({{ctrl, 'EIequation'}, fn, inputnames});

            fn = @HomogeneousBlock.setupExternalCoupling;
            inputnames = {{el, 'phi'}         , ...
                          {el, 'conductivity'}, ...
                          {ctrl, 'E'}};
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


        function state = setupEIequation(model, state)

            el   = 'ElectronicModel';
            ctrl = 'Control';

            couplingInd   = model.couplingInd;
            couplingTerms = model.couplingTerms;

            I = state.(ctrl).I;
            E = state.(ctrl).E;
            phi = state.(el).phi;

            coupterm = couplingTerms{couplingInd.input};
            faces    = coupterm.couplingfaces;
            cond_pcc = model.(el).effectiveElectronicConductivity;
            [trans_pcc, cells] = model.(el).G.getBcTrans(faces);

            state.(ctrl).EIequation = sum(cond_pcc.*trans_pcc.*(phi(cells) - E)) - I;            

        end
        
        function state = setupExternalCoupling(model, state)

            el   = 'ElectronicModel';
            ctrl = 'Control';

            E            = state.(ctrl).E;
            phi          = state.(el).phi;
            conductivity = state.(el).conductivity;
            
            couplingInd   = model.couplingInd;
            couplingTerms = model.couplingTerms;

            coupterm = couplingTerms{couplingInd.input};
            [jInput, jFaceInput] = setupExternalCoupling(model.(el), phi, E, conductivity, coupterm);

            coupterm = couplingTerms{couplingInd.output};
            [jOutput, jFaceOutput] = setupExternalCoupling(model.(el), phi, 0, conductivity, coupterm);

            state.(el).jBcSource = jInput + jOutput;

            if model.use_thermal
                state.(el).jFaceBcSource = jFaceInput + jFaceOutput;
            end
            
            state.(el).eSource = 0;

        end

        function forces = getValidDrivingForces(model)

            forces = getValidDrivingForces@PhysicalModel(model);
            forces.src = [];
            
        end


        function state = addVariables(model, state)

        % Given a state where only the primary variables are defined, this
        % functions add all the additional variables that are computed in the assembly process and have some physical
        % interpretation.
        %
        % To do so, we use getEquations function and sends dummy variable for state0, dt and drivingForces

        % Values that need to be set to get the function getEquations running

            dt = 1;
            state0 = state;
            inputparams = ControlModelInputParams([]);
            model.Control = ControlModel(inputparams);
            model.Control.controlPolicy = 'None';
            drivingForces = model.getValidDrivingForces();

            % We call getEquations to update state

            [~, state] = getEquations(model, state0, state, dt, drivingForces, 'ResOnly', true);

        end
        

        function initstate = setupInitialState(model)

            el      = 'ElectronicModel';
            thermal = 'ThermalModel';
            ctrl    = 'Control';

            nc = model.(el).G.getNumberOfCells();
            initstate.(el).phi = zeros(nc, 1);

            initstate.(ctrl).I = model.(ctrl).Imax;
            initstate.(ctrl).E = 0;
            
            if model.use_thermal
                T = model.initT;
                initstate.(thermal).T = T;
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
