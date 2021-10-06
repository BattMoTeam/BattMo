classdef BatteryInputParams 
%
% Input class for :class:`Battery <Battery.Battery>`
%
    properties
        
        
        G     % Global Grid
        SOC   % State of charge
        Ucut  % Voltage cut 
        initT % initial temperature
        
        %% parameters for the battery components
        NegativeElectrode % instance of :class:`ElectrodeInputParams`
        PositiveElectrode % instance of :class:`ElectrodeInputParams`
        Electrolyte       % instance of :class:`ElectrolyteInputParams`
        ThermalModel      % instance of :class:`ThermalModelInputParams`
                
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
