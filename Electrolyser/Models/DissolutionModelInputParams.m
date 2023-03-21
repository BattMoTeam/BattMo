classdef DissolutionModelInputParams < ComponentInputParams

    properties

        E0                     % Standard electrical potential for the dissolution reaction [V]
        c0                     % Reference concentration [mol/m^3]
        j0                     % Reference exchange current density for the dissolution reaction [A/m^2]
        MW                     % Molar mass
        rho                    % Density
        volumeFraction0        % Initial volume fraction
        volumetricSurfaceArea0 % Initial volumetric surface area
        
    end

    methods
        
        function paramobj = DissolutionModelInputParams(jsonstruct)

            paramobj = paramobj@ComponentInputParams(jsonstruct);
            
        end
        
    end
    
end
