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

        function paramobj = setup(paramobj, params)
        % params struct should contain valid fields for ComponentInputParams
        %
        % and valid fields for the methods (see implementation of those methods)
        %
        % - eac
        % - cc
        % - couplingTerms
            

            fdname = {'eac', ...
                      'cc', ...
                      'couplingTerms'};
            
            paramobj = dispatchParams(paramobj, params, 'eac');
        
            paramobj = setup@ComponentInputParams(paramobj, params);
        
        end

    end
    
end
