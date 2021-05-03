classdef ThermalComponent_ < ComponentModel
    
    methods
        
        function model = ThermalComponent_(name)
            
            model = model@ComponentModel(name);
            
            names = {'T', ...
                     'jHeatBcSource' , ...
                     'jHeatOhmSource', ...
                     'jHeatSource'   , ...
                     'jHeat'         , ...
                     'accumHeat'     , ...
                     'energyCons'};

            model.names = names;
            
            fn = @ThermalComponent.updateHeatFlux;
            inputnames = {'T'};
            fnmodel = {'.'};
            model = model.addPropFunction('jHeat', fn, inputnames, fnmodel);
            
            fn = @ThermalComponent.updateEnergyConservation;
            inputnames = {'jHeatBcSource', ...
                          'jHeatSource'  , ...
                          'jHeat'        , ...
                          'accumHeat'};
            fnmodel = {'.'};
            model = model.addPropFunction('energyCons', fn, inputnames, fnmodel);
            
        end

    end
end

