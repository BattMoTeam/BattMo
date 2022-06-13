classdef ProtonicMembraneElectrodeInputParams < ComponentInputParams
    
    properties
        
        % coefficient in Buttler-Volmer
        beta
        % Exchange current density
        iBV_0
        % Limiting current densities
        anodeLimitingCurrentDensity
        cathodeLimitingCurrentDensity
        % charge
        z
        
    end
    
    methods
        
        function paramobj = ProtonicMembraneElectrodeInputParams(jsonstruct)
            
            paramobj = paramobj@ComponentInputParams(jsonstruct);
            
        end
        
    end
    
end
