classdef ComponentInputParams 
    
    properties
        
        %% Global grid
        globG
        
        %% Cell index of current grid (indexing from global grid)
        cellind
        
        %% current Grid (contains also the mappings)
        G
                
    end
    
    methods
        
        function paramobj = ComponentInputParams(params)

            paramobj.globG = params.globG;
            paramobj.cellind = params.cellind;
            paramobj.G = genSubGrid(params.globG, params.cellind)
            
        end
        
    end
    
end
