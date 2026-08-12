classdef CO2captureChannelInputParams < ComponentInputParams

    properties

        Boundary
        
        couplingTerms

        viscosity   % viscosity
        temperature % temperature
        diameter    % tube parameter
        
        % Advanced parameter
        poiseuilleCoefficient
        
        control % Structure with fields
                % - pressure
                % - rate
                % - composition

    end

    methods
        
        function inputparams = CO2captureChannelInputParams(jsonstruct)
            
            inputparams = inputparams@ComponentInputParams(jsonstruct);

        end
        
    end

end
        
