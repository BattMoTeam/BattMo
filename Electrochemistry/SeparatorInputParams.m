classdef SeparatorInputParams < ThermalComponentInputParams

    properties
        
        thermalConductivity % (not weighted by porosity)
        thickness
        porosity
        rp
        Gurley
        
    end
    
    methods

        function paramobj = SeparatorInputParams()
            
            paramobj = paramobj@ThermalComponentInputParams();
            
        end
        
    end
    
    
end
