classdef SeaWaterElectrolyteNoPrecipitationInputParams < ElectronicComponentInputParams

    properties

        K % Reaction rates for the following reactions (depends on system)
        
        kappa % conductivity
        
        species % cell array for species 
                % Each cell is a struct with fields
                % - name (string)
                % - z (only for solutes)
                % - lambda0 (only for solutes)
                % - D : diffusion coefficients

        solutes
        solids
        logspecies
        
        quasiparticles % cell array of struct with fields
                       % - name (string)
                       % - composition : array of struct with field
                       %       - name : name of species
                       %       - coef : coefficient in the quasiParticle decomposition

        
    end
    
    methods
        
        
        function inputparams = SeaWaterElectrolyteNoPrecipitationInputParams(jsonstruct)
            
            inputparams = inputparams@ElectronicComponentInputParams(jsonstruct);
            
        end
        
    end
end
