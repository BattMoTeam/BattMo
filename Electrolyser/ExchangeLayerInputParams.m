classdef ExchangeLayerInputParams < ComponentInputParams
    
    properties
        inmrParams % struct with fields
        % - kxch : Exchange rate
        % - cT : Concentration of the relevant ion in the ionomer, mol/m^3. For a perfect membrane with only one counter-ion, this should be equal to the membrane fixed charge
        % - OH.z
        % - kML
    end
    
    methods
        
        function paramobj = ExchangeLayerInputParams(jsonstruct)

            paramobj = paramobj@ComponentInputParams(jsonstruct);
            
        end
        
    end
    
    
end
