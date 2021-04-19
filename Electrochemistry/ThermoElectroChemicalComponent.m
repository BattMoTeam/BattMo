classdef ThermoElectroChemicalComponent < ElectroChemicalComponent
    
    properties
        ThermalConductivity;
    end

    methods
        
        function model = ThermoElectroChemicalCompoenent(paramobj)
            
            model = model@ElectroChemicalComponent(paramobj);
            
            model = dispatchParams(model, paramobj, {'ThermoElectroChemicalCompoenent'});

        end

    end
end

