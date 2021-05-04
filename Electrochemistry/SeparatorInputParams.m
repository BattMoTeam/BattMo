classdef SeparatorInputParams < ComponentInputParams

    properties
        
        thickness
        porosity
        rp
        Gurley
        
        thermalConductivity
        heatCapacity
        
    end
    
    methods

        function paramobj = SeparatorInputParams()
            
            paramobj = paramobj@ComponentInputParams();
            
        end
        
    end
    
    
end
