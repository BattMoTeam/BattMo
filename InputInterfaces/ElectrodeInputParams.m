classdef ElectrodeInputParams < ComponentInputParams
    
    properties
        
        %% cell indices (in global grid)
        % shortcut used here
        % eac : ElectrodeActiveComponent
        % cc  : CurrentCollector
        eac_cellind
        cc_cellind
        
        %% coupling terms (setup by setupCouplingTerms)
        couplingTerms
        
    end
    
    methods
        
        function params = ElectrodeInputParams(params)

            params = params@ComponentInputParams(params);

            
        end
        
        function params = setupSubParams(params)
        % In this function, we setup the models ElectrodeActiveComponent and CurrentCollector
        % shortcut used here
        % eac : ElectrodeActiveComponent
        % cc  : CurrentCollector
            
            globG = params.globG;
            
            % setup ElectrodeActiveComponent
            G_eac = genSubGrid(globG, params.eac_cellind);
            params.ElectrodeActiveComponent = params.setupElectrodeActiveComponent(G_eac);
            
            % setup CurrentCollector
            G_cc = genSubGrid(globG, params.cc_cellind);
            params.CurrentCollector = CurrentCollector(G_eac);
            
        end
        
        function ElectrodeActiveComponent = setupElectrodeActiveComponent(G)
            error('virtual function')
        end
        
        function params = setupCouplingTerms(params)
        % We collect the coupling terms
            couplingTerms = {};
            couplingTerms{end + 1} = params.setupCurrentCollectorBcCoupTerm();
            couplingTerms{end + 1} = params.setupCurrentCollectorElectrodeActiveComponentCoupTerm();
            params.couplingTerms = couplingTerms;
        end

        function coupTerm = setupCurrentCollectorBcCoupTerm(params)
        % In this function, we setup coupling corresponding coupling 
            error('Virtual Function');
        end

        function coupTerm = setupCurrentCollectorElectrodeActiveComponentCoupTerm(params)
        % In this function, we setup coupling corresponding coupling 
            error('Virtual Function');
        end
        
        function params = setupModel(params)
            params.Electrode = Electrode(params);
        end

    end
    
end
