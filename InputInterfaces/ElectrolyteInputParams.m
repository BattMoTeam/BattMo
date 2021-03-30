classdef ElectrolyteInputParams < ActiveElectroChemicalComponentInputParams

    methods
        
        function paramobj = setup(paramobj, params)
            paramobj = setup@ActiveElectroChemicalComponent(paramobj, params);
        end
        
    end
    
end
