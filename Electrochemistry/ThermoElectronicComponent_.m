classdef ThermoElectronicComponent_ < ElectronicComponent_
    
    methods
        
        function model = ThermoElectronicComponent_(name)
            
            model = model@ElectronicComponent_(name);
            
            names = {model.names{:}, ...
                     'jHeatBcSource' , ...
                     'jHeatOhmSource', ...
                     'jHeatSource'   , ...
                     'jHeat'         , ...
                     'accumHeat'     , ...
                     'energyCons'};

            model.names = names;
            
            fn = @ThermoElectroChemicalComponent.updateHeatFlux;
            inputnames = {'T'};
            fnmodel = {'.'};
            model = model.addPropFunction('jHeat', fn, inputnames, fnmodel);
            
            fn = @ThermoElectroChemicalComponent.updateEnergyConservation;
            inputnames = {'jHeatBcSource', ...
                          'jHeatSource'  , ...
                          'jHeat'        , ...
                          'accumHeat'};
            fnmodel = {'.'};
            model = model.addPropFunction('energyCons', fn, inputnames, fnmodel);
            
            fn = @ThermoElectroChemicalComponent.updateOhmSource;
            inputnames = {'j'};
            fnmodel = {'.'};
            model = model.addPropFunction('jHeatOhmSource', fn, inputnames, fnmodel);
                        
            fn = @ThermoElectroChemicalComponent.updateHeatSource;
            inputnames = {'jHeatOhmSource'};
            fnmodel = {'.'};
            model = model.addPropFunction('jHeatSource', fn, inputnames, fnmodel);
            
        end

    end
end

