classdef ProtonicMembraneGasSupplyInputParams < ComponentInputParams
    
    properties
        
        molecularWeights
        permeability
        viscosity
        T
        control

        couplingTerms

    end
    
    methods
        
        function paramobj = ProtonicMembraneGasSupplyInputParams(jsonstruct)

            paramobj = paramobj@ComponentInputParams(jsonstruct);

        end
        
    end
    
end
