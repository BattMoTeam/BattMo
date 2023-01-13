classdef ExchangeLayerInputParams < ComponentInputParams
    
    properties

        kxch %
        OH % structure with field
           % - z : number of charge
        kML %
        
    end
    
    methods
        
        function paramobj = ExchangeLayerInputParams(jsonstruct)

            paramobj = paramobj@ComponentInputParams(jsonstruct);
            
        end
        
    end
    
    
end
