classdef ProtonicMembraneElectrolyteInputParams < ComponentInputParams
    
    properties

        T
        % Equilibrium H2 potential
        EH2_0
        % Equilibrium O2 potential
        EO2_0
        % Equibrium p-type conductivity 
        sigmaP_0        
        % Equibrium n-type conductivity 
        sigmaN_0
        % proton conductivity
        sigmaHp
       
    end
    
    methods
        
        function paramobj = ProtonicMembraneElectrolyteInputParams(jsonstruct)
            
            paramobj = paramobj@ComponentInputParams(jsonstruct);
            
        end
        
    end
    
end
