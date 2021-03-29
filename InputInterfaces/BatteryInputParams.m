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
        
        %% Negative electrode parameters (type ElectrodeInputParams)
        NegativeElectrodeParams        
        
        %% Positive electrode parameters (type ElectrodeInputParams)
        PositiveElectrodeParams        
                
        %% Coupling terms (describe the topological structure of the coupling)
        couplingTerms
        
        
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
        % In this function, we setup: G and the models (Electrolyte, NegativeElectrode, PositiveElectrode)
            params.NegativeElectrode = params.NegativeElectrodeParams.setupSubModels();
            params.PositiveElectrode = params.PositiveElectrodeParams.setupSubModels();
        end
        
        
        function params = setupCouplingTerms(params)
        % We collect the coupling terms
            couplingTerms = {};
            couplingTerms{end + 1} = setupNegativeElectrodeElectrolyteCoupTerm(params);
            couplingTerms{end + 1} = setupPositiveElectrodeElectrolyteCoupTerm(params);
            params.couplingTerms = couplingTerms;
            
            params.
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
