classdef MaxwellStefanGasDiffusion < MaxwellStefanDiffusion
% Model for ideal gas
    
    methods

        function model = MaxwellStefanGasDiffusion(inputparams)

            model = model@MaxwellStefanDiffusion(inputparams);
            
            model = dispatchParams(model, inputparams, fdnames);

        end
        
        function model = registerVarAndPropfuncNames(model)

            model = registerVarAndPropfuncNames@MaxwellStefanDiffusion(model);

            ncomp = model.numberOfComponents;

            fn = @MaxwellStefanGasDiffusion.updateDensity;
            inputvarnames = {VarName({}, 'massFractions', ncomp), ...
                             'pressure'};
            model = model.registerPropFunction({'density', fn, inputvarnames});


            for icomp = 1 : ncomp
                
                fn = @MaxwellStefanGasDiffusion.updateChemicalPotentials;
                outputvarnames = VarName({}, 'chemicalPotentials', ncomp, icomp);
                warning('check dependenc')
                inputvarnames  = {'temperature'                                  , ...
                                  'pressure'                                     , ...
                                  'density'                                      , ...
                                  VarName({}, 'massFractions', ncomp, icomp)};
                model = model.registerPropFunction({outputvarname, fn, inputvarnames});
                
            end
            
        end
            
        
    end
    
end

