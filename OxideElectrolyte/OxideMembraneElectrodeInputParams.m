classdef OxideMembraneElectrodeInputParams < ComponentInputParams
    
    properties
        
        T
        Eocp
        
    end
    
    methods
        
        function paramobj = OxideMembraneElectrodeInputParams(jsonstruct)
            
            paramobj = paramobj@ComponentInputParams(jsonstruct);

        end
        
    end
    
end
