classdef ThermoElectroChemicalComponent < ElectroChemicalComponent
    
    properties
        ThermalConductivity;
    end

    methods
        
        function model = ThermoElectroChemicalCompoenent(paramobj)
            
            model = model@ElectroChemicalComponent(paramobj);
            
            model = dispatchParams(model, paramobj, {'ThermoElectroChemicalCompoenent'});

        end

        function state = updateHeatFlux(model, state)
            
        end
        
        function state = updateEnergyConservation(model, state)
            
        end 
        
    end
end

