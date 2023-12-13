classdef SeaWaterBatteryInputParams < ComponentInputParams

    properties

        Cathode
        CathodeActiveMaterial
        Anode
        AnodeActiveMaterial
        Electrolyte
        
        couplingTerms    % list of coupling terms (each element is instance of couplingTerm class)

        T % Temperature is given as input (for the moment)

        include_precipitation
    end
    
    methods

        function inputparams = SeaWaterBatteryInputParams(jsonstruct)

            inputparams = inputparams@ComponentInputParams(jsonstruct);
            
            inputparams.Cathode               = HydrogenElectrodeInputParams(jsonstruct.Cathode);
            inputparams.CathodeActiveMaterial = HydrogenActiveMaterialInputParams(jsonstruct.CathodeActiveMaterial);
            inputparams.Anode                 = SeaWaterElectrodeInputParams(jsonstruct.Anode);
            inputparams.AnodeActiveMaterial   = MagnesiumActiveMaterialInputParams(jsonstruct.AnodeActiveMaterial);
            if inputparams.include_precipitation
                inputparams.Electrolyte = SeaWaterElectrolyteInputParams(jsonstruct.Electrolyte);
            else
                inputparams.Electrolyte = SeaWaterElectrolyteNoPrecipitationInputParams(jsonstruct.Electrolyte);
            end
            
            inputparams.couplingTerms = {};

            inputparams = inputparams.validateInputParams();
            
        end

        function inputparams = validateInputParams(inputparams)

            if inputparams.include_precipitation
                assert(isa(inputparams.Electrolyte, 'SeaWaterElectrolyteInputParams'), 'The input class SeaWaterElectrolyteNoPrecipitationInputParams is used. The input parameters for the electrolyte may include precipitation data.');
            end
            
        end
        
    end
    
end

