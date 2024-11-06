classdef GasDiffusionCellControlInputParams < InputParams

    properties
        
        % Cell array with the component names
        compNames

        controlElements % cell array containing struct with field
                        % - bcfaces       : index of faces which belong to same control element (indexing from GasDiffusionCell)
                        % - type          : control type (type = 1 : pressure, type = 2 : flux)
                        % - value         : value for the given control type (flux value if type = 1, pressure value if type = 2)
                        % - massFractions : mass fraction values for each component (should sum to one)
                
    end
    
    methods

        function inputparams = GasDiffusionCellControlInputParams(jsonstruct)

            inputparams = inputparams@InputParams(jsonstruct);
            
        end

        
    end
    

end
