classdef BatteryInputParams < InputParams
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
        
        function paramobj = BatteryInputParams(jsonstruct)
            
            paramobj = paramobj@InputParams(jsonstruct);
            
            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';
            elyte   = 'Electrolyte';
            thermal = 'ThermalModel';
            
            pick = @(fd) pickField(jsonstruct, fd);
            
            paramobj.(ne)      = ElectrodeInputParams(pick(ne));
            paramobj.(pe)      = ElectrodeInputParams(pick(pe));
            paramobj.(elyte)   = ElectrolyteInputParams(pick(elyte));
            paramobj.(thermal) = ThermalComponentInputParams(pick(thermal));
            paramobj.couplingTerms = {};
            
        end

    end
    
end
