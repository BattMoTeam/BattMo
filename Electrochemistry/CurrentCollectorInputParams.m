classdef CurrentCollectorInputParams < ElectronicComponentInputParams
%
% Input class for :class:`CurrentCollector <Electrochemistry.CurrentCollector>`
%
    properties
        
        couplingTerm

        thermalConductivity
        heatCapacity
        
    end
    
    methods
        
        function paramobj = CurrentCollectorInputParams();
            paramobj = paramobj@ElectronicComponentInputParams();
            paramobj.couplingTerm = struct();
            
            % we set here directly just for simplicity for the moment (hacky...)
            paramobj.EffectiveElectricalConductivity = 3.55e7;
        end
        
    end

end
