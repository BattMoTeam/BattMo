classdef CurrentCollectorInputParams < ThermoElectronicComponentInputParams

    properties
        couplingTerm;
    end
    
    methods
        
        function paramobj = CurrentCollectorInputParams();
            paramobj = paramobj@ThermoElectronicComponentInputParams();
            paramobj.couplingTerm = struct();
            
            % we set 100 here directly just for simplicity for the moment (hacky...)
            paramobj.EffectiveElectronicConductivity = 100;
        end
        
    end

end
