classdef SectorBatteryGenerator < SpiralBatteryGenerator

    
    methods
        
        function gen = SectorBatteryGenerator()
            gen = gen@SpiralBatteryGenerator();  
        end
        
        function [paramobj, gen] = setupGrid(gen, paramobj, params)
    
            gen = sectorGrid(gen);
            paramobj.G = gen.G;
            
        end

    end
    
end