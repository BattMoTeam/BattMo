classdef BatteryInputParams 
    
    properties
        
        % Global Grid
        G
        
        % Temperature and SOC
        T
        SOC
        % Input current
        J
        % Voltage cut
        Ucut
        
        
        %% parameters for the battery components
        % shortcut used here
        % ne : Negative electrode parameters (class ElectrodeInputParams)
        % pe : Positive electrode parameters (class ElectrodeInputParams)
        % elyte : Electrolyte (class ElectrolyteInputParams)
        ne;
        pe;
        elyte;        
                
        %% Coupling terms (describe the topological structure of the coupling)
        couplingTerms
        
        
    end
    
    methods
        
        function paramobj = BatteryInputParams(params)
        % params struct should contain valid fields for ComponentInputParams,
        %
        % valid fields for the methods (see implementation of those methods)
        %
        % - setupNegativeElectrode
        % - setupPositiveElectrode
        % - setupElectrolyte
        % - setupNegativeElectrodeElectrolyteCoupTerm
        % - setupPositiveElectrodeElectrolyteCoupTerm
        %
        % and fields
        %
        % - T
        % - SOC
        % - J
        % - Ucut
            
            paramobj = paramobj@setupGrid(params);
            
            paramobj = paramobj.setupVariousParams(params);
            
            % setup parameters for each of the battery component
            paramobj.ne    = paramobj.setupNegativeElectrode(params); 
            paramobj.pe    = paramobj.setupPositiveElectrode(params);
            paramobj.elyte = paramobj.setupElectrolyte(params);;

            % setup parameters for the couplings (at the battery level)
            paramobj = paramobj.setupCouplingTerms(params);
            
        end

        function paramobj = setupGrid(paramobj, params)
            error('virtual function');
        end
            
        function paramobj = setupVariousParams(paramobj, params)
            paramobj.T    = params.T;
            paramobj.SOC  = params.SOC;
            paramobj.J    = params.J;
            paramobj.Ucut = params.Ucut;
        end

        function ne_paramobj = setupNegativeElectrode(paramobj, params)
        % return object from class ElectrodeInputParams
            error('virtual function')            
        end
        
        function pe_paramobj = setupPositiveElectrode(paramobj, params)
        % return object from class ElectrodeInputParams
            error('virtual function')            
        end
        
        function elyte_paramobj = setupElectrolyte(paramobj, params)
        % return object from class ElectrolyteInputParams
            error('virtual function')            
        end
                
        function paramboj = setupCouplingTerms(paramobj, params)
        % We collect the coupling terms
            couplingTerms = {};
            couplingTerms{end + 1} = paramobj.setupNegativeElectrodeElectrolyteCoupTerm(params);
            couplingTerms{end + 1} = paramobj.setupPositiveElectrodeElectrolyteCoupTerm(params);
            paramobj.couplingTerms = couplingTerms;
        end

        function coupTerm = setupNegativeElectrodeElectrolyteCoupTerm(paramobj, params)
        % In this function, we setup coupling corresponding coupling 
            error('Virtual Function');
        end

        function coupTerm = setupPositiveElectrodeElectrolyteCoupTerm(paramobj, params)
        % In this function, we setup coupling corresponding coupling 
            error('Virtual Function');
        end

    end
    
end
