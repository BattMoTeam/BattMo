classdef ElectrolyteInputParams
    
    properties
        
        %% Global grid
        globG
        
        %% Electrode grid
        G

        %% cell indices (indexed with respect to global grid)
        cellind
        
        
    end
    
    methods
        
        function params = ElectrolyteInputParams(params)

            globG = params.globG;
            params.G = genSubGrid(globG, params.cellind)
            
        end
        
    end
    
end
