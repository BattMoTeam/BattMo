classdef CatalystLayerInputParams < ComponentInputParams
    
    properties

        j0 % Exchange current density
        inmrParams 
        % inmrParams.OH.z  : Charge
        % inmrParams.OH.c0 : OH reference concentration
        % inmrParams.E0    : standard equilibrium potential
        elyteParams
        % elyteParams.OH.c0 : reference OH concentration
        % elyteParams.E0    : standard equilibrium potential
        
        Xinmr % Fraction of specific area that is coversed with ionomer [-]
        volumetricSurfaceArea % Volumetric surface area [m^ -1]
        
    end
    
    methods
        
        function paramobj = CatalystLayerInputParams(jsonstruct)

            paramobj = paramobj@ComponentInputParams(jsonstruct);
            
        end
        
    end
    
    
end
