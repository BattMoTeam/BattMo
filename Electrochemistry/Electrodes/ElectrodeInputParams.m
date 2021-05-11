classdef ElectrodeInputParams < ComponentInputParams
    
    properties
        
        %% parameters for the electrode components
        % shortcut used here
        % eac : ElectrodeActiveComponent parameters (instance of ElectrodeActiveComponentInputParams)
        % cc  : CurrentCollector parameters (instance of CurrentCollectorInputParams)
        eac
        cc
        
        %% coupling terms
        couplingTerm
        
    end
    
    methods

        function paramobj = ElectrodeInputParams()
            paramobj = paramobj@ComponentInputParams();
            paramobj.eac = ElectrodeActiveComponentInputParams();
            paramobj.cc = CurrentCollectorInputParams();
        end

    end
    
end
