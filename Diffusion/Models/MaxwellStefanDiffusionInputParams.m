classdef MaxwellStefanDiffusionInputParams < ComponentInputParams

    properties
        
        % List of pair of diffusion coefficients. The format is a list of cells. Each cell in the list is itself a cell, where
        % - cell 1 : name of the first component (string)
        % - cell 2 : name of the second component (string)
        % - cell 3 : value of the binary diffusion coefficient between the two components
        diffusionCoefficients
        % Cell array with the component names
        componentNames
        % Struct which maps a component name (provided as a struct field name) to component index (value returned by the
        % struct)
        % Molecular weights (array with N values, where N is the number of components)
        molarWeights
        % Permeability
        permeability
        % Viscosity
        viscosity
                
    end
    
    methods

        function inputparams = MaxwellStefanDiffusionInputParams(jsonstruct)
            
            inputparams = inputparams@ComponentInputParams(jsonstruct);
            
        end
        
    end
    
end
