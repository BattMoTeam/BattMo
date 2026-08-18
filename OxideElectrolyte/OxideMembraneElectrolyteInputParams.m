classdef OxideMembraneElectrolyteInputParams < ComponentInputParams
    
    properties

        % Temperature
        T
        % Diffusion constant for hole
        Dh
        % Diffusion constant for electron
        De
        % O2- conductivity
        sigmaO2
        % Equilibrium Constant
        Keh
        % standard value of chemical electron potential
        muEl0
        
    end
    
    methods
        
        function inputparams = OxideMembraneElectrolyteInputParams(jsonstruct)
            
            inputparams = inputparams@ComponentInputParams(jsonstruct);
            
        end
        
    end
    
end
