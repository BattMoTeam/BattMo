classdef BareBatteryInputParams < InputParams
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
                
        %% Coupling terms (describe the topological structure of the coupling)
        couplingTerms
        
    end
    
    methods
        
        function paramobj = BareBatteryInputParams(jsonstruct)
            
            paramobj = paramobj@InputParams(jsonstruct);
            
            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';
            elyte   = 'Electrolyte';
            
            pick = @(fd) pickField(jsonstruct, fd);
            
            paramobj.(ne) = ElectrodeActiveComponentInputParams(pick(ne));
            paramobj.(pe) = ElectrodeActiveComponentInputParams(pick(pe));
            jsonelytestruct = pick(elyte);
            switch jsonelytestruct.electrolyteType
              case 'binary'
                paramobj.(elyte) = ElectrolyteInputParams(pick(elyte));
              otherwise
                % binary is default
                paramobj.(elyte) = ElectrolyteInputParams(pick(elyte));                
            end
        end

    end
    
end
