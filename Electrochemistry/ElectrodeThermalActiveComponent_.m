classdef ElectrodeThermalActiveComponent_ < ElectrodeActiveComponent_

    methods

        function model = ElectrodeThermalActiveComponent_(name)

            model = model@ElectroActiveChemicalComponent_(name);
            
            names = {model.names{:}, };
            model.names = names;

            
        end
        
    end
    
end

