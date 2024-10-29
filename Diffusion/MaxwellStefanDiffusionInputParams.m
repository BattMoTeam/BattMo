classdef MaxwellStefanDiffusionInputParams < InputParams

    properties
        % Diffusion maxtrix coefficients of dimension NxN (where N is the number of components)
        diffusionMatrix
        % Cell array with the component names
        compNames
        % Struct which maps a component name (provided as a struct field name) to component index (value returned by the
        % struct)
        % Molecular weights (array with N values, where N is the number of components)
        molecularWeights
        
    end
    
    methods

        function inputparams = MaxwellStefanDiffusionInputParams(jsonstruct)
            
            inputparams = inputparams@InputParams(jsonstruct);
            
        end
        
    end
    
end
