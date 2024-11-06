classdef GasDiffusionCellControlInputParams < InputParams

    properties
        
        % Cell array with the component names
        compNames

        controlElements % cell array containing struct with field
                        % - type          : control type (type = 1 : pressure, type = 2 : flux)
                        % - value         : value for the given control type (flux value if type = 1, pressure value if type = 2)
                        % - massFractions : mass fraction values for each component (should sum to one)
                        % The topology is stored in the externalCouplingTerms of the GasDiffusionCellInputParams
    end
    
    methods

        function inputparams = GasDiffusionCellControlInputParams(jsonstruct)

            inputparams = inputparams@InputParams(jsonstruct);
            
        end

        
    end
    

end
