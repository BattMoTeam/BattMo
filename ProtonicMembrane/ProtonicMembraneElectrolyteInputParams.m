classdef ProtonicMembraneElectrolyteInputParams < ComponentInputParams

    properties

        T

        sigma_n0
        sigma_prot

        dS_hyd
        dH_hyd
        Y
        Am
        dH_ox
        E_0
        steam_ratio
        Ptot
        SU
        Ea_prot
        
    end

    methods

        function paramobj = ProtonicMembraneElectrolyteInputParams(jsonstruct)

            paramobj = paramobj@ComponentInputParams(jsonstruct);

        end

    end

end
