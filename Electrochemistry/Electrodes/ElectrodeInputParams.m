classdef ElectrodeInputParams < ComponentInputParams
    
    properties
        
        %% parameters for the electrode components
        % shortcut used here
        % eac : ElectrodeActiveComponent parameters (instance of ActiveElectroChemicalComponentInputParams)
        % cc  : CurrentCollector parameters (instance of CurrentCollectorInputParams)
        eac
        cc
        
        %% coupling terms
        couplingTerm
        
    end
    
    methods

        function paramobj = ElectrodeInputParams()
            paramobj = paramobj@ComponentInputParams();
            paramobj.eac = ActiveElectroChemicalComponentInputParams();
            paramobj.cc = CurrentCollectorInputParams();
            paramobj.couplingTerm = struct();
        end

    end
    
end
