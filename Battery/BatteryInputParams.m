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

            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';
            co      = 'Coating';
            elyte   = 'Electrolyte';
            sep     = 'Separator';
            thermal = 'ThermalModel';
            ctrl    = 'Control';

            jsontruct = setDefaultJsonStructField(jsonstruct, 'include_current_collectors', false);

            jsonstruct = equalizeJsonStructFields(jsonstruct, {'include_current_collectors'      , ...
                                                               {ne, 'include_current_collectors'}, ...
                                                               {pe, 'include_current_collectors'}});

            if getJsonStructField(jsonstruct,  {pe, 'include_current_collectors'})
                jsonstruct = setDefaultJsonStructField(jsonstruct, {pe, 'use_normed_current_collector'}, true);
            end

            jsonstruct = setDefaultJsonStructField(jsonstruct, 'use_thermal', false);

            jsonstruct = equalizeJsonStructFields(jsonstruct, {{'use_thermal'}       , ...
                                                               {ne, 'use_thermal'}   , ...
                                                               {pe, 'use_thermal'}   , ...
                                                               {elyte, 'use_thermal'}, ...
                                                               {sep, 'use_thermal'}});

            jsonstruct = setDefaultJsonStructField(jsonstruct, {ne, co, 'activeMaterialModelSetup', 'composite'}, false);
            jsonstruct = setDefaultJsonStructField(jsonstruct, {pe, co, 'activeMaterialModelSetup', 'composite'}, false);

            use_thermal = getJsonStructField(jsonstruct, 'use_thermal');

            if use_thermal

                iscomp1 = getJsonStructField(jsonstruct, {ne, co, 'activeMaterialModelSetup', 'composite'});
                iscomp2 = getJsonStructField(jsonstruct, {pe, co, 'activeMaterialModelSetup', 'composite'});
                assert(~iscomp1 && ~iscomp2, 'We do not support for the moment thermal simulation for composite materials');

            end


            inputparams = inputparams@InputParams(jsonstruct);

            pick = @(fd) pickField(jsonstruct, fd);

            inputparams.(ne)      = ElectrodeInputParams(pick(ne));
            inputparams.(pe)      = ElectrodeInputParams(pick(pe));
            inputparams.(elyte)   = ElectrolyteInputParams(pick(elyte));
            inputparams.(sep)     = SeparatorInputParams(pick(sep));
            inputparams.(thermal) = ThermalComponentInputParams(pick(thermal));

            switch jsonstruct.(ctrl).controlPolicy
              case 'Impedance'
                inputparams.(ctrl) = ImpedanceControlModelInputParams(pick(ctrl));
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
              case 'timeControl'
                inputparams.(ctrl) = TimeControlModelInputParams(pick(ctrl));
              case 'TCtest'
                inputparams.(ctrl) = TCtestControlModelInputParams(pick(ctrl));
              otherwise
                error('controlPolicy %s not recognized', jsonstruct.(ctrl).controlPolicy);
            end
            inputparams.couplingTerms = {};

        end

    end

end



%{
Copyright 2021-2024 SINTEF Industry, Sustainable Energy Technology
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
