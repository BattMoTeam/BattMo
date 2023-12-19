classdef BatteryInputParams < InputParams
%
% Input parameter class for the :code:`Battery` model.
%

    properties

        G     % Computational Grid
        SOC   % Initial state of charge [-]
        initT % Initial temperature [T]

        %% parameters for the battery components

        NegativeElectrode % Negative Electrode Model, instance of :class:`Electrode <Electrochemistry.Electrodes.Electrode>`
        PositiveElectrode % Positive Electrode Model, instance of :class:`Electrode <Electrochemistry.Electrodes.Electrode>`
        Electrolyte       % Electrolyte model, instance of :class:`Electrolyte <Electrochemistry.Electrodes.Electrolyte>`
        Separator         % Separator model, instance of :class:`Separator <Electrochemistry.Electrodes.Separator>`
        ThermalModel      % Thermal model, instance of :class:`ThermalComponent <Electrochemistry.ThermalComponent>`
        Control           % Control Model

        couplingTerms % Coupling terms (describe the topological structure of the coupling between the components)

        use_thermal            % flag : true if  coupled thermal simulation should be considered
        include_current_collectors

    end

    methods

        function inputparams = BatteryInputParams(jsonstruct)

            inputparams = inputparams@InputParams(jsonstruct);

            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';
            elyte   = 'Electrolyte';
            sep     = 'Separator';
            thermal = 'ThermalModel';
            ctrl    = 'Control';

            pick = @(fd) pickField(jsonstruct, fd);

            inputparams.(ne)      = ElectrodeInputParams(pick(ne));
            inputparams.(pe)      = ElectrodeInputParams(pick(pe));
            inputparams.(elyte)   = ElectrolyteInputParams(pick(elyte));
            inputparams.(sep)     = SeparatorInputParams(pick(sep));
            inputparams.(thermal) = ThermalComponentInputParams(pick(thermal));
            
            switch jsonstruct.(ctrl).controlPolicy
              case 'CCDischarge'
                inputparams.(ctrl) = CCDischargeControlModelInputParams(pick(ctrl));
              case 'CCCharge'
                inputparams.(ctrl) = CCChargeControlModelInputParams(pick(ctrl));
              case 'CC'
                inputparams.(ctrl) = CCcontrolModelInputParams(pick(ctrl));
              case 'CCCV'
                inputparams.(ctrl) = CcCvControlModelInputParams(pick(ctrl));
              case 'powerControl'
                inputparams.(ctrl) = PowerControlModelInputParams(pick(ctrl));
              case 'CC'
                inputparams.(ctrl) = CcControlModelInputParams(pick(ctrl));
              otherwise
                error('controlPolicy %s not recognized', jsonstruct.(ctrl).controlPolicy);
            end
            inputparams.couplingTerms = {};

            inputparams = inputparams.validateInputParams();

        end

        function inputparams = validateInputParams(inputparams)

            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';
            co      = 'Coating';
            elyte   = 'Electrolyte';
            sep     = 'Separator';
            thermal = 'ThermalModel';
            ctrl    = 'Control';
            
            inputparams = mergeParameters(inputparams, {{'use_thermal'}       , ...
                                                  {ne, 'use_thermal'}   , ...
                                                  {pe, 'use_thermal'}   , ...
                                                  {elyte, 'use_thermal'}, ...
                                                  {sep, 'use_thermal'}});

            inputparams = mergeParameters(inputparams, {{'include_current_collectors'}    , ...
                                                  {ne, 'include_current_collectors'}, ...
                                                  {pe, 'include_current_collectors'}});


            inputparams.(ne)    = inputparams.(ne).validateInputParams();
            inputparams.(pe)    = inputparams.(pe).validateInputParams();
            inputparams.(elyte) = inputparams.(elyte).validateInputParams();
            inputparams.(ctrl)  = inputparams.(ctrl).validateInputParams();
            
            if inputparams.use_thermal
                inputparams.(thermal) = inputparams.(thermal).validateInputParams();

                % for the moment we do not support thermal simulation with composite material. We check for that here
                isok = strcmp(inputparams.(ne).(co).active_material_type, 'default');
                isok = isok & strcmp(inputparams.(ne).(co).active_material_type, 'default');
                assert(isok, 'We do not support for the moment thermal simulation for composite materials');

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
