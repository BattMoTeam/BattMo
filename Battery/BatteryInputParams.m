classdef BatteryInputParams < InputParams
%
% Input parameter class for the :code:`Battery` model.
%

    properties
        
        
        G     % Computational Grid
        SOC   % Initial state of charge [-]
        initT % Initial temperature [T]
        
        %% parameters for the battery components
        
        NegativeElectrode % instance of :class:`ElectrodeInputParams <Electrochemistry.ElectrodeInputParams>`
        PositiveElectrode % instance of :class:`ElectrodeInputParams <Electrochemistry.ElectrodeInputParams>`
        Electrolyte       % instance of :class:`ElectrolyteInputParams <Electrochemistry.ElectrolyteInputParams>`
        ThermalModel      % instance of :class:`ThermalComponentInputParams <Electrochemistry.ThermalComponentInputParams>`
        Control           % instance of :class:`ControlModelInputParams <Utilities.ControlModelInputParams>`
        
        couplingTerms % Coupling terms (describe the topological structure of the coupling between the components)
        
        use_thermal            % flag : true if  coupled thermal simulation should be considered
        use_particle_diffusion % flag : true if solid diffusion should be included (mainly for debugging)
        include_current_collectors

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
              case 'IEswitch'
                paramobj.(ctrl) = IEswitchControlModelInputParams(pick(ctrl));
              case 'CCCV'
                paramobj.(ctrl) = CcCvControlModelInputParams(pick(ctrl));
              case 'powerControl'
                paramobj.(ctrl) = PowerControlModelInputParams(pick(ctrl));
              case 'powerControl'
                paramobj.(ctrl) = PowerControlModelInputParams(pick(ctrl));
              case 'CC'
                paramobj.(ctrl) = CcControlModelInputParams(pick(ctrl));
              otherwise
                error('controlPolicy %s not recognized', jsonstruct.(ctrl).controlPolicy);
            end
            paramobj.couplingTerms = {};

            paramobj = paramobj.validateInputParams();
            
        end

        function paramobj = validateInputParams(paramobj)

            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';
            am      = 'ActiveMaterial';
            elyte   = 'Electrolyte';
            thermal = 'ThermalModel';
            ctrl    = 'Control';            

            paramobj = mergeParameters(paramobj, {{'use_thermal'}    , ...
                                                  {ne, 'use_thermal'}, ...
                                                  {pe, 'use_thermal'}, ...
                                                  {elyte, 'use_thermal'}});
                            
            paramobj = mergeParameters(paramobj, {{'include_current_collectors'}    , ...
                                                  {ne, 'include_current_collectors'}, ...
                                                  {pe, 'include_current_collectors'}});
                            
            paramobj = mergeParameters(paramobj, {{'use_particle_diffusion'}        , ...
                                                  {ne, am, 'use_particle_diffusion'}, ...
                                                  {pe, am, 'use_particle_diffusion'}});

            if isempty(paramobj.use_particle_diffusion)
                % set default to true
                paramobj.use_particle_diffusion = true;
                paramobj = mergeParameters(paramobj, {{'use_particle_diffusion'}        , ...
                                                      {ne, am, 'use_particle_diffusion'}, ...
                                                      {pe, am, 'use_particle_diffusion'}});
            end
            
            paramobj.(ne)    = paramobj.(ne).validateInputParams();
            paramobj.(pe)    = paramobj.(pe).validateInputParams();
            paramobj.(elyte) = paramobj.(elyte).validateInputParams();

            if paramobj.use_thermal
                paramobj.(thermal) = paramobj.(thermal).validateInputParams();
            end
            
        end

    end
    
end



%{
Copyright 2021-2023 SINTEF Industry, Sustainable Energy Technology
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
