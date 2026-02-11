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


        function state = updateDerivatives(model, state, state0, dt)

            E  = state.E;
            I  = state.I;
            E0 = state0.E;
            I0 = state0.I;

            state.dIdt = (I - I0)/dt;
            state.dEdt = (E - E0)/dt;

        end



        function state = updateControlState(model, state, state0, dt, t)

            state = updateControlState@ControlModel(model, state, state0, dt);

            datapath = fullfile(getProjectDir(), 'data', 'JP3-ckan-data', 'JP3_cell_testing_32.5-67.5pc_SoC_window_25C_SINTEF_cell_01');
            filename = 'batmax_jp3_c01_combined_QC_OCVRPT_cycling_from_SOC1.mat';
            expdata = load(fullfile(datapath, filename));
            N = 10000;
            idx = 1:1:N;
            expdata.time = expdata.time(idx);
            expdata.current = expdata.current(idx);
            state.I = interp1(expdata.time, expdata.current, t);
            fprintf('%g s: %g A (%g V)\n', t, state.I, state.E);

            state = model.updateDerivatives(state, state0, dt);

            if ~isfield(state, 'ctrlType')
                return
            end

            rsw00     = model.setupRegionSwitchFlags(state0, state0.ctrlType);
            ctrlType0 = state0.ctrlType;

            if rsw00.beforeSwitchRegion

                % We have not entered the switching region in the time step. We are not going to change control
                % in this step.
                ctrlType = ctrlType0;

            else

                % We entered the switch region in previous time step.

                currentCtrlType = state.ctrlType; % current control in the the Newton iteration
                nextCtrlType0   = model.getNextCtrlType(ctrlType0); % next control that can occur afte the previous time step control (if it changes)

                rsw0 = model.setupRegionSwitchFlags(state, ctrlType0);

                switch currentCtrlType

                  case ctrlType0

                    % The control has not changed from previous time step and we want to determine if we should change it.

                    if rsw0.afterSwitchRegion

                        % We switch to a new control because we are no longer in the acceptable region for the current
                        % control
                        ctrlType = nextCtrlType0;

                    else

                        ctrlType = ctrlType0;

                    end

                  case nextCtrlType0

                    % We do not switch back to avoid oscillation. We are anyway within the given tolerance for the
                    % control so that we keep the control as it is.

                    ctrlType = nextCtrlType0;

                  otherwise

                    error('control type not recognized');

                end

            end

            state.ctrlType = ctrlType;
            state = model.updateValueFromControl(state);

        end

        function state = updateValueFromControl(model, state)

            ImaxD = model.ImaxDischarge;
            ImaxC = model.ImaxCharge;
            Emin  = model.lowerCutoffVoltage;
            Emax  = model.upperCutoffVoltage;

            switch state.ctrlType

              case 'CC_discharge1'

                state.I = ImaxD;

              case 'CC_discharge2'

                state.I = 0;

              case 'CC_charge1'

                state.I = - ImaxC;

              case 'CV_charge2'

                state.E = Emax;

              otherwise

                error('ctrlType not recognized.')

            end

        end

        function state = updateControlEquation(model, state)

            ImaxD = model.ImaxDischarge;
            ImaxC = model.ImaxCharge;
            Emin  = model.lowerCutoffVoltage;
            Emax  = model.upperCutoffVoltage;

            E = state.E;
            I = state.I;
            ctrlType = state.ctrlType;

            switch ctrlType
              case 'CC_discharge1'
                ctrleq = I - ImaxD;
              case 'CC_discharge2'
                ctrleq = I;
              case 'CC_charge1'
                ctrleq = I + ImaxC;
              case 'CV_charge2'
                % TODO fix scaling
                ctrleq = (E - Emax)*1e5;
              otherwise
                error('ctrlType not recognized');
            end

            state.controlEquation = ctrleq;

        end

        function cleanState = addStaticVariables(model, cleanState, state)

            cleanState.numberOfCycles = state.numberOfCycles;
            cleanState.ctrlType       = state.ctrlType;

        end

        function  [arefulfilled, state] = checkConstraints(model, state, state0, dt)

            ctrlType  = state.ctrlType;
            ctrlType0 = state0.ctrlType;

            nextCtrlType = model.getNextCtrlType(ctrlType0);

            arefulfilled = true;

            rsw  = model.setupRegionSwitchFlags(state, state.ctrlType);
            rswN = model.setupRegionSwitchFlags(state, nextCtrlType);

            if (strcmp(ctrlType, ctrlType0) && rsw.afterSwitchRegion) || (strcmp(ctrlType, nextCtrlType) && ~rswN.beforeSwitchRegion)

                arefulfilled = false;

            end

        end

        function rsf = setupRegionSwitchFlags(model, state, ctrlType)

            Emin    = model.lowerCutoffVoltage;
            Emax    = model.upperCutoffVoltage;
            dIdtMin = model.dIdtLimit;
            dEdtMin = model.dEdtLimit;
            tols    = model.switchTolerances;

            E    = state.E;
            I    = state.I;

            switch ctrlType

              case 'CC_discharge1'

                before = E > Emin*(1 + tols.(ctrlType));
                after  = E < Emin*(1 - tols.(ctrlType));

              case 'CC_discharge2'

                if isfield(state, 'dEdt')
                    dEdt = state.dEdt;
                    before = abs(dEdt) > dEdtMin*(1 + tols.(ctrlType));
                    after  = abs(dEdt) < dEdtMin*(1 - tols.(ctrlType));
                else
                    before = false;
                    after  = false;
                end

              case 'CC_charge1'

                before = E < Emax*(1 - tols.(ctrlType));
                after  = E > Emax*(1 + tols.(ctrlType));

              case 'CV_charge2'

                if isfield(state, 'dIdt')
                    dIdt = state.dIdt;
                    before = abs(dIdt) > dIdtMin*(1 + tols.(ctrlType));
                    after  = abs(dIdt) < dIdtMin*(1 - tols.(ctrlType));
                else
                    before = false;
                    after  = false;
                end

              otherwise

                error('control type not recognized');

            end

            rsf = struct('beforeSwitchRegion', before, ...
                         'afterSwitchRegion' , after);

        end

        function state = updateControlAfterConvergence(model, state, state0, dt)

            initctrl = model.initialControl;

            ctrlType  = state.ctrlType;
            ctrlType0 = state0.ctrlType;

            ncycles   = state0.numberOfCycles;

            state = model.updateDerivatives(state, state0, dt);

            state.ctrlType0 = state.ctrlType;

            switch initctrl

              case 'charging'

                if ismember(ctrlType0, {'CC_discharge1', 'CC_discharge2'}) && ismember(ctrlType, {'CC_charge1', 'CV_charge2'})
                    ncycles = ncycles + 1;
                end

              case 'discharging'

                if ismember(ctrlType0, {'CC_charge1', 'CV_charge2'}) && ismember(ctrlType, {'CC_discharge1', 'CC_discharge2'})
                    ncycles = ncycles + 1;
                end

              otherwise

                error('initctrl not recognized');

            end

            state.numberOfCycles = ncycles;


        end

        function func = setupStopFunction(model)

            func = @(mainModel, state, state_prev) (state.Control.numberOfCycles >= mainModel.Control.numberOfCycles);

        end

        function func = setupControlFunction(model)
        % There is no control function in the sense that all the switches are automatically turned on and off.

            func = [];

        end

        function step = setupScheduleStep(model)

            dt = diff(model.times);
            step = struct('val', dt, ...
                          'control', ones(numel(dt), 1));

            % % Setup and a return the step structure that is part of the schedule which is used as input for
            % % :mrst:`simulateScheduleAD`. For some control type, there is a natural construction for this structure. This is
            % % why we include this method here, for convenience. It can be overloaded by derived classes. The
            % % timeSteppingParams structure by default is given by the data described in :battmofile:`Utilities/JsonSchemas/TimeStepping.schema.json`

            %     if (nargin > 1)
            %         params = timeSteppingParams;
            %     else
            %         params = [];
            %     end

            %     params = model.parseTimeSteppingStruct(params);

            %     CRate   = model.CRate;
            %     DRate   = model.DRate;
            %     ncycles = model.numberOfCycles;

            %     if isempty(ncycles)
            %         totalTime = params.totalTime;
            %     else
            %         if ~isempty(params.totalTime)
            %             warning('Both the total time and the number of cycles are given. We do not use the given total time value but compute it instead from the number of cycles.');
            %         end
            %         totalTime = ncycles*1.2*(1*hour/CRate + 1*hour/DRate);
            %     end

            %     if ~isempty(params.timeStepDuration)
            %         dt = params.timeStepDuration;
            %     else
            %         if isempty(params.numberOfTimeSteps)
            %             error('No timeStepDuration and numberOfTimeSteps are given');
            %         end
            %         n  = params.numberOfTimeSteps;
            %         dt = totalTime/n;
            %     end

            %     if params.useRampup
            %         n = params.numberOfRampupSteps;
            %     else
            %         n = 0;
            %     end

            %     dts = rampupTimesteps(totalTime, dt, n);

            %     step = struct('val', dts, 'control', ones(numel(dts), 1));

        end



        function control = setupScheduleControl(model)

            control = setupScheduleControl@ControlModel(model);
            control.CCCV2 = true;

        end

    end


    methods(Static)

        function nextCtrlType = getNextCtrlType(ctrlType)

            switch ctrlType

              case 'CC_discharge1'

                nextCtrlType = 'CC_discharge2';

              case 'CC_discharge2'

                nextCtrlType = 'CC_charge1';

              case 'CC_charge1'

                nextCtrlType = 'CV_charge2';

              case 'CV_charge2'

                nextCtrlType = 'CC_discharge1';

              otherwise

                error('ctrlType not recognized.')

            end

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
