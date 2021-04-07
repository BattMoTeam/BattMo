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

        
        function paramobj = setupGrid(paramobj, params)
            
            globG = params.globG;
            cellind = params.cellind;
            
            paramobj.globG = globG;
            paramobj.cellind = cellind;
            paramobj.G = genSubGrid(globG, cellind);
        
        end
        
        
    end
    
end
