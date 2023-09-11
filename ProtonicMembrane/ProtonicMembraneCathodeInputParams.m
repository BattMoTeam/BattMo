classdef ProtonicMembraneCathodeInputParams < ProtonicMembraneElectrodeInputParams
    
    properties

        R
        
    end
    
    methods
        
        function paramobj = ProtonicMembraneCathodeInputParams(jsonstruct)
            
            paramobj = paramobj@ProtonicMembraneElectrodeInputParams(jsonstruct);
            
        end
        
    end
    
end
