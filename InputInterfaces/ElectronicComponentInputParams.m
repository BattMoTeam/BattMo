classdef ElectronicComponentInputParams < ComponentInputParams

    properties
        EffectiveElectronicConductivity;
    end
    
    methods
        
        function paramobj = ElectronicComponentInputParams(params)
        % params struct should contain valid fields for ComponentInputParams
        %
        % and the field
        %
        % - EffectiveElectronicConductivity
            paramobj = paramobj@ComponentInputParams(params);
            paramobj.EffectiveElectronicConductivity = params.EffectiveElectronicConductivity;
        end
        
    end
    
end
