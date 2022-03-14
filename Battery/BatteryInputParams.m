classdef BatteryInputParams < InputParams
%
% Input parameter class for :class:`Battery <Battery.Battery>`
%

    properties
        
        
        G     % Computational Grid
        SOC   % Initial state of charge [-]
        initT % Initial temperature [T]
        
        %% parameters for the battery components
        NegativeElectrode % instance of :class:`ElectrodeInputParams <Electrochemistry.ElectrodeInputParams>`
        PositiveElectrode % instance of :class:`ElectrodeInputParams <Electrochemistry.ElectrodeInputParams>`
        Electrolyte       % instance of :class:`ElectrolyteInputParams <Electrochemistry.ElectrolyteInputParams>`
        ThermalModel      % instance of :class:`ThermalModelInputParams <Electrochemistry.ThermalModelInputParams>`
        Control           % instance of :class:`ControlModelInputParams <Electrochemistry.ControlModelInputParams>`    
        
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
            ctrl    = 'Control';            

            pick = @(fd) pickField(jsonstruct, fd);
            
            paramobj.(ne)      = ElectrodeInputParams(pick(ne));
            paramobj.(pe)      = ElectrodeInputParams(pick(pe));
            paramobj.(elyte)   = ElectrolyteInputParams(pick(elyte));
            paramobj.(thermal) = ThermalComponentInputParams(pick(thermal));
            switch jsonstruct.(ctrl).controlPolicy
              case 'EIswitch'
                paramobj.(ctrl)    = ControlModelInputParams(pick(ctrl));
              case 'CCCV'
                paramobj.(ctrl)    = CcCvControlModelInputParams(pick(ctrl));
            end
            paramobj.couplingTerms = {};
            
        end

    end
    
end



%{
Copyright 2021-2022 SINTEF Industry, Sustainable Energy Technology
and SINTEF Digital, Mathematics & Cybernetics.

This file is part of The Battery Modeling Toolbox BattMo

BattMo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BattMo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BattMo.  If not, see <http://www.gnu.org/licenses/>.
%}
