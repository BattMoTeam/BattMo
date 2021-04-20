classdef BatteryInputParams 
    
    properties
        
        % Global Grid
        G
        %  SOC
        SOC
        % Input current
        J
        % Voltage cut
        Ucut
        % initial temperature
        initT
        
        %% parameters for the battery components
        % shortcut used here
        % ne : Negative electrode parameters (instance of ElectrodeInputParams)
        % pe : Positive electrode parameters (instance of ElectrodeInputParams)
        % elyte : Electrolyte (instance ElectrolyteInputParams)
        ne;
        pe;
        elyte;        
                
        %% Coupling terms (describe the topological structure of the coupling)
        couplingTerms
        
    end
    
    methods
        
        function paramobj = BatteryInputParams()
            paramobj.ne = ElectrodeInputParams();
            paramobj.pe = ElectrodeInputParams();
            paramobj.elyte = ElectrolyteInputParams();
            parmobj.couplingTerms = {};
        end

    end
    
end
