classdef IonomerMembraneInputParams < ElectronicComponentInputParams
    
    properties
        
        volumeFraction
        
        H2O % with fields
        % H2O.c0 : Reference concentration
        % H2O.D : diffusion coefficient for water
        % H2O.MW : Molar mass (needed for function groupHydration which is only needed in setup of initial condition and not for assembly)

        
        OH % with fields
        % OH.xi : OH occupation
        % OH.z : Charge number
        % OH.t : Transference number

        cT % Total concentration of charged group (one scalar value)
        
        V % molar volume (needed for function groupHydration which is only needed in setup of initial condition and not for assembly)

        tortuosity % cell-valued coefficient
        
    end
    
    methods
        
        function paramobj = IonomerMembraneInputParams(jsonstruct)
            
            paramobj = paramobj@ElectronicComponentInputParams(jsonstruct);
            
        end
        
    end
    
    
end
