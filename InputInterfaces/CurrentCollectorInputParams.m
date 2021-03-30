classdef CurrentCollectorInputParams < ElectronicComponentInputParams
% It is now an instantiation (we set the EffectiveElectronicConductivity here). Done for simplicity.
    methods
        
        function paramobj = CurrentCollectorInputParams();
            paramobj = paramobj@ElectronicComponentInputParams();
            paramobj.EffectiveElectronicConductivity = 100;
        end
        
    end

end
