classdef ProtonicMembraneElectrodeInputParams < ComponentInputParams
    
    properties
        
        T
        E_0
        
    end
    
    methods
        
        function inputparams = ProtonicMembraneElectrodeInputParams(jsonstruct)
            
            inputparams = inputparams@ComponentInputParams(jsonstruct);
            
        end
        
    end
    
end
