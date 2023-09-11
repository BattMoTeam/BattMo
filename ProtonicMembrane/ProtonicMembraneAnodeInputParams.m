classdef ProtonicMembraneAnodeInputParams < ProtonicMembraneElectrodeInputParams
    
    properties
        
        % coefficient in Buttler-Volmer
        beta
        
        % charge-transfer current densigy
        i0

        % Limiting current densities
        ila % Anode
        ilc % Cathode

        Rct
        ila0
        n
        SU
        pH2O_in
        
    end
    
    methods
        
        function paramobj = ProtonicMembraneAnodeInputParams(jsonstruct)
            
            paramobj = paramobj@ProtonicMembraneElectrodeInputParams(jsonstruct);
            
        end
        
    end
    
end
