classdef ElectrodeInputParams < ComponentInputParams
%
% Input class for :class:`Electrode <Electrochemistry.Electrodes.Electrode>`
%        
    properties
        
        %% parameters for the electrode components
        ElectrodeActiveComponent
        CurrentCollector
        
        %% coupling terms
        couplingTerm
        
    end
    
    methods

        function paramobj = ElectrodeInputParams(jsonstruct)

            paramobj = paramobj@ComponentInputParams(jsonstruct);

            eac = 'ElectrodeActiveComponent';
            cc  = 'CurrentCollector';
            
            pick = @(fd) pickField(jsonstruct, fd);
            paramobj.(eac) = ElectrodeActiveComponentInputParams(pick(eac));
            paramobj.(cc) = CurrentCollectorInputParams(pick(cc));
            
        end

    end
    
end
