classdef GasDiffusionCell < MaxwellStefanGasDiffusion

    properties

        Control
        
    end
    
    methods

        function model = GasDiffusionCell(inputparams)

            model = model@MaxwellStefanGasDiffusion(inputparams);

            model.Control = GasDiffusionCellControl(inputparams.Control);

            model.subModelNameList{end + 1} = 'Control';
            
        end

        function model = registerVarAndPropfuncNames(model)

            ncomp = model.numberOfComponents;
            
            model = registerVarAndPropfuncNames@MaxwellStefanGasDiffusion(model);

            model = model.setAsStaticVarName(VarName({}, 'sources', ncomp));
            
            fn = @GasDiffusionCell.updateControlValue;
            model = model.registerPropFunction({{'Control', 'value'}, fn, {{'Control', 'type'}}});

            fn = @GasDiffusionCell.updateControlFluxBoundaryEquation;
            inputvarnames = {VarName({'Boundary'}, 'diffusionFluxes', ncomp), ...
                             VarName({'Control'}, 'flux')};
            outputvarnames = VarName({'Control'}, 'fluxBoundaryEquations', ncomp);
            model = model.registerPropFunction({outputvarnames, fn, inputvarnames});

            fn = @GasDiffusionCell.updateBoundaryMassFractions;
            inputvarnames = {VarName({'Control'}, 'massFractions', ncomp)};
            outputvarnames = VarName({'Boundary'}, 'massFractions', ncomp);
            model = model.registerPropFunction({outputvarnames, fn, inputvarnames});
            
            fn = @GasDiffusionCell.updateControlBoundaryPressureEquation;
            inputvarnames = {{'Control', 'pressure'}, ...
                             {'Boundary', 'pressure'}};
            outputvarname = {'Control', 'pressureBoundaryEquation'};
            model = model.registerPropFunction({outputvarname, fn, inputvarnames});

        end
        
    end

end
