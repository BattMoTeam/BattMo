classdef CO2captureInputParams < InputParams
    
    properties

        Feed
        Permeate
        permeances % structure with permeability
        length
        
    end

    methods

        function inputparams = CO2captureInputParams(jsonstruct)

            jsonstruct = equalizeStructFields(jsonstruct, {'length'      , ...
                                                           {'Feed', 'length'}, ...
                                                           {'Permeate', 'length'}});

            jsonstruct = equalizeStructFields(jsonstruct, {'gasSpecies'      , ...
                                                           {'Feed', 'gasSpecies'}, ...
                                                           {'Permeate', 'gasSpecies'}});
            
            inputparams = inputparams@InputParams(jsonstruct);

            inputparams.Feed     = CO2captureChannelInputParams(jsonstruct.Feed);
            inputparams.Permeate = CO2captureChannelInputParams(jsonstruct.Permeate);
            
        end

    end

end
