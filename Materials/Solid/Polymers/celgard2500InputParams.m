classdef celgard2500InputParams < SeparatorInputParams
    
    methods
        
        function paramobj = celgard2500InputParams()
            
            paramobj.thickness = 10e-6;
            paramobj.porosity  = 0.55;
            paramobj.rp        = 0.064e-6 ./ 2;
            paramobj.Gurley    = 200;

        end
        
    end
    
end

