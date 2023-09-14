classdef OxideMembraneElectrodeInputParams < ComponentInputParams
    
    properties
        
        T
        Rct   % charge transfer resistance
        muEl0 % standard value of chemical electron potential
        pO2   % O2 pressure used to compute Eocp
        
    end
    
    methods
        
        function paramobj = OxideMembraneElectrodeInputParams(jsonstruct)
            
            paramobj = paramobj@ComponentInputParams(jsonstruct);

        end
        
    end
    
end
