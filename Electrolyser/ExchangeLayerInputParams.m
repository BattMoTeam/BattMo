classdef ExchangeLayerInputParams < ComponentInputParams
    
    properties
        inmrParams % struct with fields
        % kxch
        % cT
        % OH.z
        % kML        
    end
    
    methods
        
        function paramobj = ExchangeLayerInputParams(jsonstruct)

            paramobj = paramobj@ComponentInputParams(jsonstruct);
            
        end
        
    end
    
    
end
