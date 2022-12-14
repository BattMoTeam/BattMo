classdef IonomerMembraneInputParams < ElectronicComponentInputParams
    
    properties
        
        liquidVolumeFraction
        
        H2O % with fields
        % H2O.c0 :  Reference concentration
        % H2O.D :  diffusion coefficient for water
        
        OH % with fields
        % OH.xi : OH occupation
        % OH.z :  
        % OH.t :        
        
    end
    
    methods
        
        function paramobj = IonomerMembraneInputParams(jsonstruct)
            
            paramobj = paramobj@ElectronicComponentInputParams(jsonstruct);
            
        end
        
    end
    
    
end
