classdef ProtonicMembraneElectrodeInputParams < ComponentInputParams
    
    properties
        
        T
        E_0
        
    end
    
    methods
        
        function paramobj = ProtonicMembraneElectrodeInputParams(jsonstruct)
            
            paramobj = paramobj@ComponentInputParams(jsonstruct);
            
        end
        
    end
    
end
