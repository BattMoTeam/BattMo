classdef ThermoElectroChemicalComponent_ < ElectroChemicalComponent_
    
    methods
        
        function model = ThermoElectroChemicalCompoenent_(name)
            
            model = model@ElectroChemicalComponent_(name);
            
            names = {model.names{:}, ...
                     'jHeatBCSource', ...
                     'jHeatSource'  , ...
                     'jHeat'        , ...
                     'energyCons'};

            model.names = names;
            
            fn = @ThermoElectroChemicalComponent.updateHeatFlux;
            inputnames = {'T'};
            fnmodel = {'.'};
            model = model.addPropFunction('jHeat', fn, inputnames, fnmodel);
            
            fn = @ThermoElectroChemicalComponent.updateEnergyConservation;
            inputnames = {'jHeatBCSource', ...
                          'jHeatSource'  , ...
                          'jHeat'};
            fnmodel = {'.'};
            
            model = model.addPropFunction('energyCons', fn, inputnames, fnmodel);
            

        end

    end
end

