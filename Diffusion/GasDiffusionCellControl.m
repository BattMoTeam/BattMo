classdef GasDiffusionCellControl < BaseModel

    properties
        
        % Cell array with the component names
        compNames
        % Struct which maps a component name (provided as a struct field name) to component index (value returned by the
        % struct)
        compInds
        % Number of components
        numberOfComponents

        controlValues
        
    end
    
    methods

        function model = GasDiffusionCellControl(inputparams)

            model = model@BaseModel();
            
            fdnames = {'compNames'};

            model = dispatchParams(model, inputparams, fdnames);
            
            model.numberOfComponents = numel(model.compNames);

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
            % Total flux : one for each control boundary
            varnames{end + 1} = 'flux';
            % Pressure : one value for each control boundary
            varnames{end + 1} = 'pressure';
            % Mass Fraction : one value for each control boundary
            varnames{end + 1} = VarName({}, 'massFractions', ncomp);
            % Control value : one value for each control boundary
            varnames{end + 1} = 'value';            
            % Control flux equations : one value per component, for each control boundary
            varnames{end + 1} = VarName({}, 'fluxBoundaryEquations', ncomp);
            % Control pressure equation : one valuefor each control boundary
            varnames{end + 1} = 'pressureBoundaryEquation';
            % Control types : one value for each control boundary. Two types
            % - pressure control : type = 1
            % - flux control : type = 1
            varnames{end + 1} = 'type';
            % Control value equations : one value for each control boundary
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
