classdef FlatBatteryGenerator < SpiralBatteryGenerator

    properties
        depth
    end
    
    methods
        
        function gen = FlatBatteryGenerator()
            gen = gen@SpiralBatteryGenerator();  
        end
        
        function [paramobj, gen] = setupGrid(gen, paramobj, params)
    
            gen = flatGrid(gen);
            paramobj.G = gen.G;
            
        end

        function paramobj = updateBatteryInputParams(gen, paramobj, params)
                    
            gen.nwindings = params.nwindings;
            gen.depth     = params.depth;
            gen.widthDict = params.widthDict;
            gen.nrDict    = params.nrDict;
            gen.nas       = params.nas;
            gen.L         = params.L;
            gen.nL        = params.nL;
            
            paramobj = gen.setupBatteryInputParams(paramobj, []);
            
        end
        
    end
    
end