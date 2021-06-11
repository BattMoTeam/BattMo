classdef BatteryInputParams 
    
    properties
        
        % Global Grid
        G
        %  SOC
        SOC
        % Voltage cut
        Ucut
        % initial temperature
        initT
        
        %% parameters for the battery components
        NegativeElectrode % (instance of ElectrodeInputParams)
        PositiveElectrode % (instance of ElectrodeInputParams)
        Electrolyte       % (instance ElectrolyteInputParams)
        ThermalModel      % (instance ThermalModelInputParams)
                
        %% Coupling terms (describe the topological structure of the coupling)
        couplingTerms
        
    end
    
    methods
        
        function paramobj = BatteryInputParams()
            
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            elyte = 'Electrolyte';
            thermal = 'ThermalModel';
            
            paramobj.(ne) = ElectrodeInputParams();
            paramobj.(pe) = ElectrodeInputParams();
            paramobj.(elyte) = ElectrolyteInputParams();
            paramobj.(thermal) = ThermalComponentInputParams();
            paramobj.couplingTerms = {};
            
        end

    end
    
end
