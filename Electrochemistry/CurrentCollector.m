classdef CurrentCollector < ElectronicComponent
% Current collector model 
    
    properties
        volumeFraction
        electronicConductivity
    end
    
    methods
        
        function model = CurrentCollector(params)
            
            model = model@ElectronicComponent(params);
            
            model.electronicConductivity = params.electronicConductivity;
            model.volumeFraction = params.volumeFraction;
            
        end
    end
    
end
                             
