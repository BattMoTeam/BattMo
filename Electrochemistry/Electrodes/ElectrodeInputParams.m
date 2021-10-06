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

        function paramobj = ElectrodeInputParams()

            eac = 'ElectrodeActiveComponent';
            cc  = 'CurrentCollector';

            paramobj = paramobj@ComponentInputParams();
            
            paramobj.(eac) = ElectrodeActiveComponentInputParams();
            paramobj.(cc) = CurrentCollectorInputParams();
            
        end

    end
    
end
