classdef ProtonicMembraneGasSupplyInputParams < ComponentInputParams
    
    properties
        
        molecularWeights
        permeability
        viscosity
        diffusionCoefficients % vector of diffusion coefficients, one per component
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
