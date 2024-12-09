classdef ProtonicMembraneCathodeInputParams < ProtonicMembraneElectrodeInputParams
    
    properties

        R_ct_H
        Ptot
        
    end
    
    methods
        
        function inputparams = ProtonicMembraneCathodeInputParams(jsonstruct)
            
            inputparams = inputparams@ProtonicMembraneElectrodeInputParams(jsonstruct);
            
        end
        
    end
    
end
