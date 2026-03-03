classdef OxideMembraneElectrodeInputParams < ComponentInputParams
    
    properties
        
        T
        Rct   % charge transfer resistance
        muEl0 % standard value of chemical electron potential
        pO2   % O2 pressure used to compute Eocp
        Keh   % Equilibrium constant for the hole-electron reaction
        
    end
    
    methods
        
        function inputparams = OxideMembraneElectrodeInputParams(jsonstruct)
            
            inputparams = inputparams@ComponentInputParams(jsonstruct);

        end
        
    end
    
end
