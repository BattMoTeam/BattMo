classdef ElectronicComponentInputParams < ComponentInputParams

    properties
        EffectiveElectronicConductivity;
    end
    
    methods
        

        function paramobj = setup(paramobj, params)
        % params struct should contain valid fields for ComponentInputParams
        %
        % and the field
        %
        % - EffectiveElectronicConductivity
            paramobj = setup@ComponentInputParams(paramobj, params);
            paramobj.EffectiveElectronicConductivity = params.EffectiveElectronicConductivity;        
        end
        
    
    end
    
end
