classdef CcCvControlModel2 < ControlModel
%
% Constant-Current-Constant-Voltage Control
%


    properties

        times
        values

        %
        % Charge and discharge rates
        CRate
        DRate

        lowerCutoffVoltage
        upperCutoffVoltage

        dIdtLimit
        dEdtLimit

        numberOfCycles

        % Control used initially. String that can take one of the following values
        % - 'discharging'
        % - 'charging'
        initialControl

        % This values are initiated depending on C/D rate values and battery model
        ImaxCharge
        ImaxDischarge

        %
        switchTolerances

    end


    methods

        function model = CcCvControlModel2(inputparams)

            model = model@ControlModel(inputparams);

            fdnames = {'CRate'             , ...
                       'DRate'             , ...
                       'lowerCutoffVoltage', ...
                       'upperCutoffVoltage', ...
                       'dEdtLimit'         , ...
                       'dIdtLimit'         , ...
                       'numberOfCycles'    , ...
                       'initialControl'    , ...
                       'switchTolerances'  , ...
                       'times'             , ...
                       'values'};

            model = dispatchParams(model, inputparams, fdnames);

            if isempty(model.numberOfCycles)
                warning('Number of cycles has not been given in CCCV control. We use numberOfCycles = 1.');
                model.numberOfCycles = 1;
            end

        end

        function model = registerVarAndPropfuncNames(model)

            model = registerVarAndPropfuncNames@ControlModel(model);

            varnames = {};
            % Control type : string that can take following value
            % - CC_discharge1
            % - CC_discharge2
            % - CC_charge1
            % - CV_charge2
            varnames{end + 1} = 'ctrlType';

            % Estimate of the time derivative for voltage and current
            varnames{end + 1} = 'dIdt';
            varnames{end + 1} = 'dEdt';

            % Variable to store number of cycles
            varnames{end + 1} = 'numberOfCycles';
            model = model.registerVarNames(varnames);

            % The following variables are not used in the residual assembly

            varnames = {'numberOfCycles', ...
                        'dIdt'          , ...
                        'dEdt'};

            model = model.setAsStaticVarNames(varnames);
            model = model.setAsExtraVarNames(varnames);

            % Register the functions

            fn = @CcCvControlModel2.updateControlEquation;
            model = model.registerPropFunction({'controlEquation', fn, {'ctrlType', 'E', 'I'}});

            fn = @CcCvControlModel2.updateDerivatives;
            fn = {fn, @(prop) PropFunction.accumFuncCallSetupFn(prop)};
            inputvarnames = {'E', 'I'};
            model = model.registerPropFunction({'dIdt', fn, {'E', 'I'}});
            model = model.registerPropFunction({'dEdt', fn, {'E', 'I'}});

        end

        function ctrlVal = computeInput(model, t)

            ctrlVal = interp1(model.times, model.values, t);

        end


        function state = updateControlState(model, state, state0, dt)

            state = updateControlState@ControlModel(model, state, state0, dt);

            % state = model.updateDerivatives(state, state0, dt);

            Emin = model.lowerCutoffVoltage;
            Emax = model.upperCutoffVoltage;

            %ctrlType = 'current';
            ctrlType = state.ctrlType;

            if state.E <= Emin && strcmp(state.ctrlType, 'current')
                ctrlType = 'Emin';
                % state.E = Emin;
                %keyboard;
            end

            if state.E > Emax && strcmp(state.ctrlType, 'current')
                ctrlType = 'Emax';
                state.E = Emax;
            end

            if strcmp(state.ctrlType, 'Emin') && state.I > state.ctrlVal
                ctrlType = 'current';
                state.I = state.ctrlVal;
            end

            if strcmp(state.ctrlType, 'Emax') && state.I < state.ctrlVal
                ctrlType = 'current';
                state.I = state.ctrlVal;
            end

            state.ctrlType = ctrlType;


        end

        function state = updateControlEquation(model, state)

            Emin  = model.lowerCutoffVoltage;
            Emax  = model.upperCutoffVoltage;

            E = state.E;
            I = state.I;
            ctrlType = state.ctrlType;

            switch ctrlType
              case 'current'
                ctrleq = I - state.ctrlVal;
              case 'Emin'
                ctrleq = (E - Emin)*1e5;
              case 'Emax'
                ctrleq = (E - Emax)*1e5;
              otherwise
                error('ctrlType %s not recognized', ctrlType);
            end

            state.controlEquation = ctrleq;

        end

        function cleanState = addStaticVariables(model, cleanState, state)

            %cleanState.numberOfCycles = state.numberOfCycles;
            cleanState.ctrlType       = state.ctrlType;

        end

        function  [arefulfilled, state] = checkConstraints(model, state, state0, dt)

            % ctrlType  = state.ctrlType;
            % ctrlType0 = state0.ctrlType;

            % nextCtrlType = model.getNextCtrlType(ctrlType0);

            % arefulfilled = true;

            % rsw  = model.setupRegionSwitchFlags(state, state.ctrlType);
            % rswN = model.setupRegionSwitchFlags(state, nextCtrlType);

            % if (strcmp(ctrlType, ctrlType0) && rsw.afterSwitchRegion) || (strcmp(ctrlType, nextCtrlType) && ~rswN.beforeSwitchRegion)

            %     arefulfilled = false;

            % end

            arefulfilled = true;

        end

        function state = updateControlAfterConvergence(model, state, state0, dt)

            initctrl = model.initialControl;

            ctrlType  = state.ctrlType;
            ctrlType0 = state0.ctrlType;

            %ncycles   = state0.numberOfCycles;

            state.ctrlType0 = state.ctrlType;

            %state.numberOfCycles = ncycles;


        end

        function func = setupControlFunction(model)
        % There is no control function in the sense that all the switches are automatically turned on and off.

            func = [];

        end

        function step = setupScheduleStep(model)

            dt = diff(model.times);
            step = struct('val', dt, 'control', ones(numel(dt), 1));

        end



        function control = setupScheduleControl(model)

            control = setupScheduleControl@ControlModel(model);
            control.CCCV2 = true;

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
