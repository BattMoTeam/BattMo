classdef ThermoElectroChemicalComponent_ < ElectroChemicalComponent_
    
    methods
        
        function model = ThermoElectroChemicalCompoenent_(name)
            
            model = model@ElectroChemicalComponent_(name);
            
            names = {model.names{:}, ...
                    }

        end

    end
end

