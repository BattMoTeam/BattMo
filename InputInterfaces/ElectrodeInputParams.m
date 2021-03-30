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


        function paramobj = setup(paramobj, params)
        % params struct should contain valid fields for ComponentInputParams
        %
        % and valid fields for the methods (see implementation of those methods)
        %
        % - eac
        % - cc
        % - couplingTerms
            

            paramobj.eac = getparam(params, 'eac');
            paramobj.cc  = getparam(params, 'cc');
            paramobj.couplingTerms = getparam(params, 'couplingTerms');
        
            paramobj = setup@ComponentInputParams(paramobj, params);
        
        end

    end
    
end
