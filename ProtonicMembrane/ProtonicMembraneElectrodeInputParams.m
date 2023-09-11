classdef ProtonicMembraneElectrodeInputParams < ComponentInputParams
    
    properties
        
        T
        Eocp
        
    end
    
    methods
        
        function paramobj = ProtonicMembraneElectrodeInputParams(jsonstruct)
            
            paramobj = paramobj@ComponentInputParams(jsonstruct);
            
        end
        
    end
    
end
