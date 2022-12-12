classdef PorousTransportLayerInputParams < ElectronicComponentInputParams
    
    properties
        
        solidVolumeFraction % Solid volume fraction
        leverettCoefficient % Leverett coefficient that enters in the computation of the capillary pressure
        theta % Water contact angle, enters in the computation of the capillary pressure
        permeability % Permeability [Darcy]
        BruggemanCoefficient % Bruggeman coefficient
        sp % species struct 
        % sp.OH.MW    : Molecular weight [kg mol^-1]
        % sp.OH.V0    : Partial molar volume
        % sp.OH.D     : Diffustion coefficient
        % sp.OH.t     : Transference coefficient 
        % sp.OH.z     : Charge number
        % sp.K.MW
        % sp.K.V0
        % sp.H2O.MW   : 
        % sp.H2O.beta : Water equilibrium constant
        % sp.H2O.kLV  : liquid-vapor exchange rate
        % sp.H2O.mu0  : Standard chemical potential
        
        
    end
    
    methods
        
        function paramobj = HydrogenPorousTransportLayerInputParams();
            paramobj = paramobj@ElectronicComponentInputParams();
            
        end
        
    end
    
    
end
