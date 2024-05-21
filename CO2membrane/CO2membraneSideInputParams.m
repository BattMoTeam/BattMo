classdef CO2membraneSideInputParams < ComponentInputParams

    properties

        Boundary
        
        couplingTerm

        viscosity   % viscosity
        temperature % temperature
        diameter    % tube parameter
        
        % Advanced parameter
        poiseuilleCoefficient
        
    end

    methods
        
        function inputparams = CO2membraneSide(jsonstruct)
            
            inputparams = inputparams@ComponentInputParams(jsonstruct);

            inputparams.Boundary = CO2membraneSideBoundaryInputParams(jsonstruct.Boundary);


        end
        
    end

end
        
