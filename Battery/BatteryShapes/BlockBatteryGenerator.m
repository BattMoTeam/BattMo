classdef BlockBatteryGenerator < SpiralBatteryGenerator
    
    methods
        
        function gen = BlockBatteryGenerator()
            gen = gen@SpiralBatteryGenerator();  
        end
        

        function [paramobj, gen] = setupGrid(gen, paramobj, params)
    
            gen = flatGrid(gen);
            paramobj.G = gen.G;
            
        end

    end
    
end
