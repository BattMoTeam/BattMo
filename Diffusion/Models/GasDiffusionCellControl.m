classdef GasDiffusionCellControl < BaseModel

    properties
        
        % Cell array with the component names
        componentNames
        % Struct which maps a component name (provided as a struct field name) to component index (value returned by the
        % struct)
        compInds
        % Number of components
        numberOfComponents

        controlElements % cell array containing struct with field
                        % - bcfaces       : index of faces which belong to same control element (indexing from GasDiffusionCell)
                        % - type          : control type (type = 1 : pressure, type = 2 : flux)
                        % - value         : value for the given control type (flux value if type = 1, pressure value if type = 2)
                        % - massFractions : mass fraction values for each component (should sum to one)
        
        nControls % number of control elements (equal to size of controlElements);

        mappings % structure with following fields
                 % - ctrlToBc : mapping from control values to boundary values (control values are vector with nControls values)
                 % - bcToCtrl : mapping from boundary values to control values. It is the transpose of the ctrlToBc and
                 %              will sum up the boundary values that belong to the same control element
                 % - fluxToValue     : maps the flux vector variables to the values that corresponds to the value variable
                 %                     (see method updateValueEquation and the definition of the variables in
                 %                     registerVarAndPropfuncNames)
                 % - pressureToValue : same as fluxToValue but for the pressure vector variable
    end
    
    methods

        function model = GasDiffusionCellControl(inputparams)

            model = model@BaseModel();
            
            fdnames = {'componentNames', ...
                       'controlElements'};

            model = dispatchParams(model, inputparams, fdnames);
            
            model.numberOfComponents = numel(model.componentNames);
            model.nControls          = numel(model.controlElements);

            model = model.setupValueMappings();
            
            ncomp = model.numberOfComponents;
            
            for icomp = 1 : ncomp

                name = model.componentNames{icomp};
                compInds.(name) = icomp;
                
            end

            model.compInds = compInds;
            
        end


        function model = setupValueMappings(model)
        % We setup the mappings fluxToValue and pressureToValue

            typetbl.type = cellfun(@(elt) elt.type, model.controlElements);
            typetbl = IndexArray(typetbl);

            fluxtbl.type = 2;
            fluxtbl = IndexArray(fluxtbl);
            
            pressuretbl.type = 1;
            pressuretbl = IndexArray(pressuretbl);

            map = TensorMap();
            map.fromTbl  = fluxtbl;
            map.toTbl    = typetbl;
            map.mergefds = {'type'};
            map = map.setup();

            fluxToValue = map.getMatrix();

            map = TensorMap();
            map.fromTbl  = pressuretbl;
            map.toTbl    = typetbl;
            map.mergefds = {'type'};
            map = map.setup();

            pressureToValue = map.getMatrix();
            
            model.mappings.fluxToValue     = fluxToValue;
            model.mappings.pressureToValue = pressureToValue;
            
        end

        function model = registerVarAndPropfuncNames(model)

            model = registerVarAndPropfuncNames@BaseModel(model);

            ncomp = model.numberOfComponents;

            varnames = {};
            % Total flux : one for each control boundary element
            varnames{end + 1} = 'flux';
            % Pressure : one value for each control boundary element
            varnames{end + 1} = 'pressure';
            % Mass Fraction : one value for each control boundary element
            varnames{end + 1} = VarName({}, 'massFractions', ncomp);
            % Control value : one value for each control boundary element
            varnames{end + 1} = 'value';            
            % Control flux equations : one value per component, for each control boundary element
            varnames{end + 1} = VarName({}, 'fluxBoundaryEquations', ncomp);
            % Control pressure equation : one valuefor each control boundary element
            varnames{end + 1} = 'pressureBoundaryEquation';
            % Control types : one value for each control boundary element. Two types
            % - pressure control : type = 1
            % - flux control : type = 1
            varnames{end + 1} = 'type';
            % Control value equations : one value for each control boundary element
            % The flux or pressure variables are set to the value, depending on the type
            varnames{end + 1} = 'valueEquation';

            model = model.registerVarNames(varnames);

            model = model.setAsStaticVarName('type');

            fn = @GasDiffusionCellControl.updateValueEquation;
            inputvarnames = {'value', 'flux', 'pressure'};
            model = model.registerPropFunction({'valueEquation', fn, inputvarnames});
            
        end

        function state = updateValueEquation(model, state)

            m = model.mappings;
            
            state.valueEquation =  state.value - m.fluxToValue*state.flux - m.pressureToValue*state.pressure;
            
        end
        
    end
    

end
