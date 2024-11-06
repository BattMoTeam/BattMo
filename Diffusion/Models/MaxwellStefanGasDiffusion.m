classdef MaxwellStefanGasDiffusion < MaxwellStefanDiffusion
% Model for ideal gas
    
    methods

        function model = MaxwellStefanGasDiffusion(inputparams)

            model = model@MaxwellStefanDiffusion(inputparams);

        end
        
        function model = registerVarAndPropfuncNames(model)

            model = registerVarAndPropfuncNames@MaxwellStefanDiffusion(model);

            ncomp = model.numberOfComponents;

            fn = @MaxwellStefanGasDiffusion.updateDensity;
            inputvarnames = {VarName({}, 'massFractions', ncomp), ...
                             'pressure'};
            model = model.registerPropFunction({'density', fn, inputvarnames});

            fn = @MaxwellStefanGasDiffusion.updateDensity;
            inputvarnames = {VarName({'Boundary'}, 'massFractions', ncomp), ...
                             VarName({'Boundary'}, 'pressure')};
            outputvarname = VarName({'Boundary'}, 'density');
            model = model.registerPropFunction({outputvarname, fn, inputvarnames});

            for icomp = 1 : ncomp
                
                fn = @MaxwellStefanGasDiffusion.updateChemicalPotentials;
                outputvarname = VarName({}, 'chemicalPotentials', ncomp, icomp);
                % warning('check dependency')
                inputvarnames  = {'temperature', ...
                                  'pressure'   , ...
                                  'density'    , ...
                                  VarName({}, 'massFractions', ncomp, icomp)};
                model = model.registerPropFunction({outputvarname, fn, inputvarnames});

                fn = @MaxwellStefanGasDiffusion.updateChemicalPotentials;
                outputvarname = VarName({'Boundary'}, 'chemicalPotentials', ncomp, icomp);
                % warning('check dependency')
                inputvarnames  = {VarName({'Boundary'}, 'temperature'), ...
                                  VarName({'Boundary'}, 'pressure')   , ...
                                  VarName({'Boundary'}, 'density')    , ...
                                  VarName({'Boundary'}, 'massFractions', ncomp, icomp)};
                model = model.registerPropFunction({outputvarname, fn, inputvarnames});
                                
            end
            
        end
            
        
    end
    
end

