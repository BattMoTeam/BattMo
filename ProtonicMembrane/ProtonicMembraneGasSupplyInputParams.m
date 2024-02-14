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
        
        function inputparams = ProtonicMembraneGasSupplyInputParams(jsonstruct)

            inputparams = inputparams@ComponentInputParams(jsonstruct);

        end
        
    end
    
end
