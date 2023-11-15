classdef ProtonicMembraneAnodeInputParams < ProtonicMembraneElectrodeInputParams
    
    properties
        
        % coefficient in Buttler-Volmer
        beta
        
        % Limiting current densities
        ila % Anode
        ilc % Cathode

        n
        
        R_ct_0
        Ea_ct      
        SU         
        O2_conc_feed
        steam_ratio
        Ptot

        gasSupplyType
        
    end
    
    methods
        
        function paramobj = ProtonicMembraneAnodeInputParams(jsonstruct)
            
            paramobj = paramobj@ProtonicMembraneElectrodeInputParams(jsonstruct);
            
        end
        
    end
    
end
