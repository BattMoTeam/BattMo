classdef BatteryInputParams
    
    properties
        
        % Temperature and SOC
        T
        SOC

        % Input current
        J
        
        % Voltage cut
        Ucut
        
        % Grid of the whole model (parent grid)
        G

        %% Electrolyte model
        Electrolyte
        
        %% Negative electrode model
        NegativeElectrode
        
        %% Positive electrode model
        PositiveElectrode
        
        %% Separator model
        sep
        
        %% Current collector for negative electrode model
        NegativeCurrentCollector
        
        %% Current collector for positive electrode model 
        PositiveCurrentCollector
        
        %% Coupling terms (describe the topological structure of the coupling)

        % electrolyte - negative electrode (location where chemical reactions take place)
        coupTermNegativeElectrodeElectrolyte
        
        % electrolyte - positive electrode (location where chemical reactions take place)
        coupTermPositiveElectrodeElectrolyte 

        % current collector and electrode (current continuity is imposed at the faces)
        coupTermNegativeCurrentCollectorNegativeElectrode  
        coupTermPositiveCurrentCollectorPositiveElectrode  

        % Boundary conditions for the current collectors
        coupTermNegativeCurrentCollectorBc  
        coupTermPositiveCurrentCollectorBc  
        
    end
    
    methods
        
        function params = BatteryInputParams(params)

            params = params.setupVariousParams();
            params = params.setupSubModels();
            params = params.setupCouplingTerms();
            
        end
        
        function params = setupVariousParams(param)
        % In this function, we define :  SOC, T,  J and  Ucut
            error('virtual function')
        end
        
        function params = setupSubModels(params)
        % In this function, we setup: G and the models (Electrolyte, sep, NegativeElectrode, PositiveElectrode, NegativeCurrentCollector, PositiveCurrentCollector)
            error('Virtual Function');
        end
        
        function params = setupCouplingTerms(params)
            params.coupTermNegativeElectrodeElectrolyte = setupNegativeElectrodeElectrolyteCoupTerm(params);
            params.coupTermPositiveElectrodeElectrolyte = setupPositiveElectrodeElectrolyteCoupTerm(params);
            params.coupTermNegativeCurrentCollectorNegativeElectrode  = setupNegativeCurrentCollectorNegativeElectrodeCoupTerm(params);
            params.coupTermPositiveCurrentCollectorPositiveElectrode  = setupPositiveCurrentCollectorPositiveElectrodeCoupTerm(params);
            params.coupTermNegativeCurrentCollectorBc  = setupNegativeCurrentCollectorBcCoupTerm(params);
            params.coupTermPositiveCurrentCollectorBc  = setupPositiveCurrentCollectorBcCoupTerm(params);
        end

        function coupTerm = setupNegativeCurrentCollectorBcCoupTerm(params)
        % In this function, we setup coupling corresponding coupling 
            error('Virtual Function');
        end

        function coupTerm = setupNegativeCurrentCollectorNegativeElectrodeCoupTerm(params)
        % In this function, we setup coupling corresponding coupling 
            error('Virtual Function');
        end
        
        function coupTerm = setupPositiveCurrentCollectorBcCoupTerm(params)
        % In this function, we setup coupling corresponding coupling 
            error('Virtual Function');
        end

        function coupTerm = setupPositiveCurrentCollectorPositiveElectrodeCoupTerm(params)
        % In this function, we setup coupling corresponding coupling 
            error('Virtual Function');            
        end

        function coupTerm = setupNegativeElectrodeElectrolyteCoupTerm(params)
        % In this function, we setup coupling corresponding coupling 
            error('Virtual Function');
        end

        function coupTerm = setupPositiveElectrodeElectrolyteCoupTerm(params)
        % In this function, we setup coupling corresponding coupling 
            error('Virtual Function');
        end

    end
    
end
