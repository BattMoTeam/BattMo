classdef ProtonicMembraneCathodeInputParams < ProtonicMembraneElectrodeInputParams
    
    properties

        R_ct_H
        Ptot
        
    end
    
    methods
        
        function paramobj = ProtonicMembraneCathodeInputParams(jsonstruct)
            
            paramobj = paramobj@ProtonicMembraneElectrodeInputParams(jsonstruct);
            
        end
        
    end
    
end
