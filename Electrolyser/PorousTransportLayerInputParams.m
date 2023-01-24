classdef PorousTransportLayerInputParams < ElectronicComponentInputParams
    
    properties
        
        solidVolumeFraction  % Solid volume fraction [-]
        leverettCoefficients  % Leverett coefficient that enters in the computation of the capillary pressure, see leverett.m 
        theta                % Water contact angle, enters in the computation of the capillary pressure [-]
        permeability         % Permeability [Darcy]
        BruggemanCoefficient % Bruggeman coefficient [-]
        sp % species struct 
        % sp.OH.MW    : Molecular weight of OH [kg mol^-1]
        % sp.OH.V0    : Partial molar volume of OH [m^3 mol^-1]
        % sp.OH.D     : Diffusion coefficient [m^2 s^-1]
        % sp.OH.t     : Transference coefficient [-]
        % sp.OH.z     : Charge number [-]
        % sp.K.MW     : Molecular weight of K [kg mol^-1]
        % sp.K.V0     : Partial molar volume of K [m^3 mol^-1]
        % sp.H2O.MW   : Molecular weight of H2O [kg mol^-1]
        % sp.H2O.kLV  : Liquid-vapor exchange rate
        % sp.H2O.mu0  : Standard chemical potential
        % sp.H2O.V0   : Partial molar volume of H2O [m^3 mol^-1]
        externalCouplingTerm

        Boundary % structure with field
                 % sp : the structure contains the molecular weight of the gas components.
        
    end
    
    methods
        
        function paramobj = PorousTransportLayerInputParams(jsonstruct)
            
            paramobj = paramobj@ElectronicComponentInputParams(jsonstruct);
            
        end
        
    end
    
    
end
