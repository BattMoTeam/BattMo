classdef CatalystLayerInputParams < ComponentInputParams
    
    properties

        j0    % Exchange current density
        E0    % Standard equilibrium potential
        Eref  % Reference potential

        sp % species struct with field
        % - OH.z  : Charge
        % - OH.c0 : OH reference concentration

        n % Number of electron transfer
        
        alpha                 % coefficient in the exponent in Butler-Volmer equation [-]
        Xinmr                 % Fraction of specific area that is coversed with ionomer [-]
        volumetricSurfaceArea % Volumetric surface area [m^ -1]
        
    end
    
    methods
        
        function paramobj = CatalystLayerInputParams(jsonstruct)

            paramobj = paramobj@ComponentInputParams(jsonstruct);
            
        end
        
    end
    
    
end
