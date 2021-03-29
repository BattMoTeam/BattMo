classdef ElectronicComponentInputParams < ComponentInputParams

    properties
        EffectiveElectronicConductivity;
    end
    
    methods
        
        function paramobj = ElectronicComponentInputParams(params)
            paramobj = params@ComponentInputParams(params);
            paramobj.EffectiveElectronicConductivity = params.EffectiveElectronicConductivity;
        end
        
    end
    
end
