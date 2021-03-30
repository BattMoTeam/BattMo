classdef ElectrodeInputParams < ComponentInputParams
    
    properties
        
        %% parameters for the electrode components
        % shortcut used here
        % eac : ElectrodeActiveComponent parameters (class ActiveElectroChemicalComponentInputParams)
        % cc  : CurrentCollector parameters (class CurrentCollectorInputParams)
        eac
        cc
        
        %% coupling terms
        couplingTerms
        
    end
    
    methods

        function paramobj = ElectrodeInputParams()
            paramobj.eac = ActiveElectroChemicalComponentInputParams();
            paramobj.cc = CurrentCollectorInputParams();
            paramobj.couplingTerms = {};
        end

    end
    
end
