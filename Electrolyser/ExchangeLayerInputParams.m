classdef ExchangeLayerInputParams < ComponentInputParams
    
    properties
        
        kxch % Exchange rate
        OH   % structure with field
             % - z : number of charge
        kML  % Ionomer sorption coefficient
        
    end
    
    methods
        
        function paramobj = ExchangeLayerInputParams(jsonstruct)

            paramobj = paramobj@ComponentInputParams(jsonstruct);
            
        end
        
    end
    
    
end
