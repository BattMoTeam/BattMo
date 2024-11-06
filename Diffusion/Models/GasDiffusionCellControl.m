classdef GasDiffusionCellControl < BaseModel

    properties
        
        % Cell array with the component names
        compNames
        % Struct which maps a component name (provided as a struct field name) to component index (value returned by the
        % struct)
        compInds
        % Number of components
        numberOfComponents

        controlElements % cell array containing struct with field
                        % - bcfaces       : index of faces which belong to same control element (indexing from GasDiffusionCell)
                        % - type          : control type (type = 1 : pressure, type = 2 : flux)
                        % - values        : value for the given control type (flux value if type = 1, pressure value if type = 2)
                        % - massFractions : mass fraction values for each component (should sum to one)
        
        nControls % number of control elements (equal to size of controlElements);
        
    end
    
    methods

        function model = GasDiffusionCellControl(inputparams)

            model = model@BaseModel();
            
            fdnames = {'compNames', ...
                       'controlElements'};

            model = dispatchParams(model, inputparams, fdnames);
            
            model.numberOfComponents = numel(model.compNames);
            model.nControls          = numel(model.controlElements);
            
            ncomp = model.numberOfComponents;
            
            for icomp = 1 : ncomp

                name = model.compNames{icomp};
                compInds.(name) = icomp;
                
            end

            model.compInds = compInds;
            
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
        
        
    end
    

end
