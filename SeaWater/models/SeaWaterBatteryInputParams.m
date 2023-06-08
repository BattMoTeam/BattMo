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

        function paramobj = SeaWaterBatteryInputParams(jsonstruct)

            paramobj = paramobj@ComponentInputParams(jsonstruct);
            
            paramobj.Cathode               = HydrogenElectrodeInputParams(jsonstruct.Cathode);
            paramobj.CathodeActiveMaterial = HydrogenActiveMaterialInputParams(jsonstruct.CathodeActiveMaterial);
            paramobj.Anode                 = SeaWaterElectrodeInputParams(jsonstruct.Anode);
            paramobj.AnodeActiveMaterial   = MagnesiumActiveMaterialInputParams(jsonstruct.AnodeActiveMaterial);
            if paramobj.include_precipitation
                paramobj.Electrolyte = SeaWaterElectrolyteInputParams(jsonstruct.Electrolyte);
            else
                paramobj.Electrolyte = SeaWaterElectrolyteNoPrecipitationInputParams(jsonstruct.Electrolyte);
            end
            
            paramobj.couplingTerms = {};

            paramobj = paramobj.validateInputParams();
            
        end

        function paramobj = validateInputParams(paramobj)

            if paramobj.include_precipitation
                assert(isa(paramobj.Electrolyte, 'SeaWaterElectrolyteInputParams'), 'The input class SeaWaterElectrolyteNoPrecipitationInputParams is used. The input parameters for the electrolyte may include precipitation data.');
            end
            
        end
        
    end
    
end

