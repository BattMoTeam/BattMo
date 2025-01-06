classdef ThermalModel < BaseModel

    methods

        function model = registerVarAndPropfuncNames(model)
            
            model = registerVarAndPropfuncNames@BaseModel(model);
            
            varnames = {};
            % Temperature
            varnames{end + 1} = 'T';
            % Accumulation term
            varnames{end + 1} = 'accumTerm';
            % Flux flux
            varnames{end + 1} = 'flux';
            % Heat source
            varnames{end + 1} = 'source';
            % Energy conservation
            varnames{end + 1} = 'energyCons';
            
            model = model.registerVarNames(varnames);
            
            fn = @ThermalModel.updateFlux;
            model = model.registerPropFunction({'flux', fn, {'T'}});

            fn = @ThermalModel.updateAccumTerm;
            model = model.registerPropFunction({'accumTerm', fn, {'T'}});
            
            fn = @ThermalModel.updateEnergyCons;
            inputnames = {'accumTerm', 'flux', 'source'};            
            model = model.registerPropFunction({'energyCons', fn, inputnames});

            
        end

    end
    
end
