classdef ProtonicMembraneElectrodeInputParams < ComponentInputParams
    
    properties
        
        T
        
    end
    
    methods
        
        function paramobj = ProtonicMembraneElectrodeInputParams(jsonstruct)
            
            paramobj = paramobj@ComponentInputParams(jsonstruct);
            
        end
        
    end
    
end
