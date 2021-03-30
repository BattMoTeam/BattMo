classdef CurrentCollector < ElectronicComponent
% Current collector model : At the moment equivalent to an ElectronicComponent
    
    methods
        
        function model = CurrentCollector(paramobj)
            
            model = model@ElectronicComponent(paramobj);
            
            % The parameter EffectiveElectronicConductivity in CurrentCollectorInputParams is given as scalar
            
            model.EffectiveElectronicConductivity = model.EffectiveElectronicConductivity*ones(model.G.cells.num, 1);
            
        end
        
    end
    
end
                             
