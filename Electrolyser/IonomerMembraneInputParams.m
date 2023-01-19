classdef IonomerMembraneInputParams < ElectronicComponentInputParams
    
    properties
        
        volumeFraction
        
        H2O % with fields
        % H2O.c0 : Reference concentration
        % H2O.D : diffusion coefficient for water
        
        OH % with fields
        % OH.xi : OH occupation
        % OH.z : Charge number
        % OH.t : Transference number

        cT % Total concentration of charged group
        
    end
    
    methods
        
        function paramobj = IonomerMembraneInputParams(jsonstruct)
            
            paramobj = paramobj@ElectronicComponentInputParams(jsonstruct);
            
        end
        
    end
    
    
end
