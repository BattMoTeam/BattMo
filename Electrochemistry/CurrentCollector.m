classdef CurrentCollector < ElectrochemicalComponent
% Current collector model 
    
    properties
        volumeFraction
        electronicConductivity
    end
    
    methods
        function model = CurrentCollector(G, cells)
            
            model = model@ElectrochemicalComponent();
            
            % setup grid
            G = genSubGrid(G, cells);
            model.G = G;
            
            % setup discrete differential operators
            model.operators = localSetupOperators(G);

            model.electronicConductivity = 100;
            model.volumeFraction = ones(G.cells.num, 1);
            model.EffectiveElectronicConductivity = (model.electronicConductivity) .* (model.volumeFraction).^1.5;            

            
        end
    end
    
end
                             
