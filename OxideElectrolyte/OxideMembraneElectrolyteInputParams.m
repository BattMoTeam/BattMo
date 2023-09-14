classdef OxideMembraneElectrolyteInputParams < ComponentInputParams
    
    properties

        T
        % Equilibrium H2 potential
        EH2_0
        % Equilibrium O2 potential
        EO2_0

        sigmaN_0
        Y
        dH_hyd
        dS_hyd
        Ea_prot
        pH2O_in
        pH2O_neg
        SU
        t_p_O2

    end
    
    methods
        
        function paramobj = OxideMembraneElectrolyteInputParams(jsonstruct)
            
            paramobj = paramobj@ComponentInputParams(jsonstruct);
            
        end
        
    end
    
end
