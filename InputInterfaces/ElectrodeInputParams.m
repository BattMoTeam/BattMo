classdef ElectrodeInputParams < ComponentInputParams
    
    properties
        
        %% parameters for the electrode components
        % shortcut used here
        % eac : ElectrodeActiveComponent parameters (class ActiveElectroChemicalComponentInputParams)
        % cc  : CurrentCollector parameters (class CurrentCollectorInputParams)
        eac
        cc
        
        %% coupling terms (setup by setupCouplingTerms)
        couplingTerms
        
    end
    
    methods
        
        function paramobj = ElectrodeInputParams(params)
        % params struct should contain valid fields for ComponentInputParams
        %
        % and valid fields for the methods (see implementation of those methods)
        %
        % - setupElectrodeActiveComponent
        % - setupCurrentCollector
        % - setupCurrentCollectorBcCoupTerm
        % - setupCurrentCollectorElectrodeActiveComponentCoupTerm
            
            paramobj = paramobj@ComponentInputParams(params);
            paramobj.eac = paramobj.setupElectrodeActiveComponent(params);
            paramobj.cc  = paramobj.setupCurrentCollector(params);
            paramobj.couplingTerms = paramobj.setupCouplingTerms(params);
            
        end
        
        function paramobj = setupCouplingTerms(paramobj, params)
        % We collect the coupling terms
            couplingTerms = {};
            couplingTerms{end + 1} = paramobj.setupCurrentCollectorBcCoupTerm(params);
            couplingTerms{end + 1} = paramobj.setupCurrentCollectorElectrodeActiveComponentCoupTerm(params);
            paramobj.couplingTerms = couplingTerms;
        end

        function eac_paramobj = setupElectrodeActiveComponent(paramobj, params)
        % return object from class ElectrodeActiveComponentInputParams
            error('virtual function')            
        end
        
        function eac_paramobj = setupCurrentCollector(paramobj, params)
        % return object from class CurrentCollectorInputParams
            error('virtual function')            
        end
        
        function coupTerm = setupCurrentCollectorBcCoupTerm(paramobj, params)
        % In this function, we setup coupling corresponding coupling 
            error('Virtual Function');
        end

        function coupTerm = setupCurrentCollectorElectrodeActiveComponentCoupTerm(paramobj, params)
        % In this function, we setup coupling corresponding coupling 
            error('Virtual Function');
        end

    end
    
end
