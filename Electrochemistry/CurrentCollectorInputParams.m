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
        
        function paramobj = CurrentCollectorInputParams(jsonstruct);
            paramobj = paramobj@ElectronicComponentInputParams(jsonstruct);
            paramobj.couplingTerm = struct();
        end
        
    end

end
