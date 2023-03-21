classdef CatalystLayerInputParams < ComponentInputParams
    
    properties

        j0    % Exchange current density
        E0    % Standard equilibrium potential
        Eref  % Reference potential

        sp % species struct with field
        % - OH.z  : Charge
        % - OH.c0 : OH reference concentration

        n % Number of electron transfer
        
        alpha                  % coefficient in the exponent in Butler-Volmer equation [-]
        Xinmr                  % Fraction of specific area that is coversed with ionomer [-]
        volumetricSurfaceArea0 % Volumetric surface area [m^ -1]

        tortuosity % Tortuosity [-]

        include_dissolution % True if dissolution model is included (if not given, it is set to false)

        DissolutionModel
        
    end
    
    methods
        
        function paramobj = CatalystLayerInputParams(jsonstruct)

            paramobj = paramobj@ComponentInputParams(jsonstruct);

            if isempty(paramobj.include_dissolution)
                paramobj.include_dissolution = false;
            elseif paramobj.include_dissolution
                paramobj.DissolutionModel = DissolutionModelInputParams(jsonstruct.DissolutionModel);
            end

            paramobj = paramobj.validateInputParams();
            
        end

        function paramobj = validateInputParams(paramobj)

            paramobj = validateInputParams@ComponentInputParams(paramobj);
            
            dm = 'DissolutionModel';

            if paramobj.include_dissolution
                paramobj = mergeParameters(paramobj, {{'volumetricSurfaceArea0'}, {dm, 'volumetricSurfaceArea0'}});
                paramobj.(dm) = paramobj.(dm).validateInputParams();
            end
            

        end

    end
    
    
end
