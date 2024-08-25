classdef CcCvRestControlModel < ControlModel
%
% Constant-Current-Constant-Voltage Control
%

    
    properties

        %
        % Charge and discharge rates
        CRate
        DRate
        
        lowerCutoffVoltage
        upperCutoffVoltage

        Ilimit

        totalRestTime
        
        % This values are initiated depending on C/D rate values and battery model
        ImaxCharge
        ImaxDischarge
        
        %
        tolerances
        
    end
    
    
    methods

        function model = CcCvRestControlModel(inputparams)

            model = model@ControlModel(inputparams);
            
            fdnames = {'CRate'             , ...
                       'DRate'             , ...
                       'lowerCutoffVoltage', ...
                       'upperCutoffVoltage', ...
                       'Ilimit'            , ...
                       'totalRestTime'};
            
            model = dispatchParams(model, inputparams, fdnames);

            % values of these relative tolerances should be smaller than 1
            tolerances = struct('CC_charge1'   , 1e-3, ...
                                'CV_charge2'   , 1e-3, ...
                                'CC_discharge1', 1e-3, ...
                                'CV_discharge2', 1e-3, ...
                                'Rest'         , 1e-3);
            
            model.tolerances = tolerances;
            
        end

        function model = registerVarAndPropfuncNames(model)

            model = registerVarAndPropfuncNames@ControlModel(model);
            
            varnames = {};
            % Control type : string that can take following value
            % - CC_charge1
            % - CV_charge2
            % - Rest
            % - CC_discharge1
            % - CV_discharge2
            varnames{end + 1} = 'ctrlType';
            varnames{end + 1} = 'restTime';

            model = model.registerVarNames(varnames);

            model = model.setAsStaticVarName('restTime');
            
            % Register the functions
            fn = @CcCvRestControlModel.updateControlEquation;
            model = model.registerPropFunction({'controlEquation', fn, {'ctrlType', 'E', 'I'}});

        end

        function state = updateRestTime(model, state, state0, dt)

            switch state0.ctrlType
                
              case 'CC_charge1'
                
                state.restTime = 0;
                
              case 'CV_charge2'
                
                if strcmp(state.ctrlType, 'Rest')
                    state.restTime = dt;
                else
                    state.restTime = 0;
                end
                
              case 'Rest'
                
                if strcmp(state.ctrlType, 'Rest')
                    state.restTime = state0.restTime + dt;
                else
                    % do nothing
                end
                
              case {'CC_discharge1', 'CV_discharge2'}
                
                % do nothing
                
              otherwise
                error('ctrlType not recognized');
                
            end

        end

        
        function state = updateControlState(model, state, state0, dt)
        % Called after each Newton iteration

            state = updateControlState@ControlModel(model, state, state0, dt);

            % A new state has been computed with a given control type (stored in state.ctrlType). In this function we find
            % the new control type that is the best candidate. We use an overall conservative approach
            
            currentCtrlType = state.ctrlType;
            ctrlType0       = state0.ctrlType; 

            nextCtrlType0 = model.getNextCtrlType(ctrlType0);
            
            rsw0 = model.setupRegionSwitchFlags(state, ctrlType0);
            
            switch currentCtrlType

              case nextCtrlType0
                
                % We are in the situation where the control has switched: The current control is the next one after the
                % previous control type. We want to determine if we switch back. 
                
                if  rsw0.beforeSwitchRegion
                    
                    % We switch back because we are before the beginning of the switching region
                    ctrlType = ctrlType0;

                else

                    % We keep the control type as it is
                    ctrlType = nextCtrlType0;

                end

              case ctrlType0

                % The control has not changed from previous and we want to determine if we should change it. 

                if rsw0.afterSwitchRegion

                    % We switch to a new control because we are no longer in the acceptable region for the current
                    % control
                    ctrlType = nextCtrlType0;
                    
                else

                    ctrlType = ctrlType0;

                end

              otherwise

                error('control type not recognized');

            end

            state.ctrlType = ctrlType;
            
            state = model.updateValueFromControl(state);

        end

        function state = updateValueFromControl(model, state)

            ImaxC = model.ImaxCharge;
            ImaxD = model.ImaxDischarge;
            Emax  = model.upperCutoffVoltage;
            Emin  = model.lowerCutoffVoltage;
            
            switch state.ctrlType
                
              case 'CC_charge1'

                state.I = - ImaxC;
                
              case 'CV_charge2'

                state.E = Emax;

              case 'Rest'

                state.I = 0;

              case 'CC_discharge1'

                state.I = ImaxD;

              case 'CV_discharge2'

                state.E = Emin;

              otherwise
                
                error('ctrlType not recognized.')

            end
            
        end
        
        function state = updateControlEquation(model, state)
            
            ImaxC = model.ImaxCharge;
            ImaxD = model.ImaxDischarge;
            Emin  = model.lowerCutoffVoltage;
            Emax  = model.upperCutoffVoltage;
            
            E = state.E;
            I = state.I;            
            ctrlType = state.ctrlType;
            
            switch ctrlType
              case 'CC_charge1'
                ctrleq = I + ImaxC;
              case 'CV_charge2'
                % TODO fix scaling
                ctrleq = (E - Emax)*1e5;
              case 'Rest'
                ctrleq = I;
              case 'CC_discharge1'
                ctrleq = I - ImaxD;
              case 'CV_discharge2'
                ctrleq = (E - Emin)*1e5;
              otherwise
                error('ctrlType not recognized');
            end
            
            state.controlEquation = ctrleq;
            
        end

        function cleanState = addStaticVariables(model, cleanState, state)

            cleanState.ctrlType = state.ctrlType;
            cleanState.restTime = state.restTime;
            
        end
        
        function  [arefulfilled, state] = checkConstraints(model, state, state0, dt)

            ctrlType  = state.ctrlType;
            ctrlType0 = state0.ctrlType;
            
            nextCtrlType = model.getNextCtrlType(ctrlType0);

            arefulfilled = true;
            
            rsw = model.setupRegionSwitchFlags(state, state.ctrlType);

            if strcmp(ctrlType, ctrlType0) && rsw.afterSwitchRegion
                
                arefulfilled = false;
                state.ctrlType = nextCtrlType;
                state = model.updateValueFromControl(state);

            end
                
        end
        
        function rsf = setupRegionSwitchFlags(model, state, ctrlType)

            Emin       = model.lowerCutoffVoltage;
            Emax       = model.upperCutoffVoltage;
            Ilimit     = model.Ilimit; % positive value
            totaltRest = model.totalRestTime;
            
            tols    = model.tolerances;
            
            E     = state.E;
            I     = state.I;
            tRest = state.restTime;

            switch ctrlType
                
              case 'CC_charge1'
                
                before = (E - Emax)/Emax < -tols.(ctrlType);
                after  = (E - Emax)/Emax > tols.(ctrlType);

              case 'CV_charge2'

                % Recall that in charge phase, I is negative by convention
                before = (-I - Ilimit)/Ilimit > tols.(ctrlType);
                after  = (-I - Ilimit)/Ilimit < -tols.(ctrlType);

              case 'Rest'

                before = (tRest - totaltRest)/totaltRest < -tols.(ctrlType);
                after  = (tRest - totaltRest)/totaltRest > tols.(ctrlType);
                
              case 'CC_discharge1'

                before = (E - Emin)/Emin > tols.(ctrlType);
                after  = (E - Emin)/Emin < -tols.(ctrlType);

              case 'CV_discharge2'
                                
                before = (E - Emin)/Emin > tols.(ctrlType);
                after  = (E - Emin)/Emin < -tols.(ctrlType);
                
              otherwise

                error('control type not recognized');
            end

            rsf = struct('beforeSwitchRegion', before, ...
                         'afterSwitchRegion' , after);
            
        end

        function state = updateControlAfterConvergence(model, state, state0, dt)

            state = updateControlAfterConvergence@ControlModel(model, state, state0, dt);
            state = model.updateRestTime(state, state0, dt);
            state = model.updateControlState(state, state0, dt);
            
        end        
        
        function func = setupStopFunction(model)

            function dostop = stopfunction(mainModel, state, state_prev)

                dostop = false;

                if strcmp(state.Control.ctrlType, 'CV_discharge2') 
                    dostop = true;
                end
                
            end
            
            func = @stopfunction;
            
        end

        function func = setupControlFunction(model)
        % There is no control function in the sense that all the switches are automatically turned on and off.
            
            func = [];
            
        end

        function step = setupScheduleStep(model, timeSteppingParams)
            
            if (nargin > 1)
                params = timeSteppingParams;
            else
                params = [];
            end

            params = model.parseTimeSteppingStruct(params);
            
            CRate   = model.CRate;
            DRate   = model.DRate;
            trest   = model.totalRestTime;
            
            if isempty(params.totalTime) 
                totalTime = 1.5*(1*hour/CRate + 1*hour/DRate + trest);
            end

            if ~isempty(params.timeStepDuration)
                dt = params.timeStepDuration;
            else
                if isempty(params.numberOfTimeSteps) 
                    error('No timeStepDuration and numberOfTimeSteps are given');
                end
                n  = params.numberOfTimeSteps;
                dt = totalTime/n;
            end

            if params.useRampup
                n = params.numberOfRampupSteps;
            else
                n = 0;
            end

            dts = rampupTimesteps(totalTime, dt, n);
            
            step = struct('val', dts, 'control', ones(numel(dts), 1));

        end
        
        function control = setupScheduleControl(model)
            
            control = setupScheduleControl@ControlModel(model);
            control.CCCVrest = true;
            
        end
        
    end
    
    methods(Static)

        function nextCtrlType = getNextCtrlType(ctrlType)

            switch ctrlType
                
              case 'CC_charge1'
                
                nextCtrlType = 'CV_charge2';

              case 'CV_charge2'
                
                nextCtrlType = 'Rest';

              case 'Rest'

                nextCtrlType = 'CC_discharge1';
                
              case 'CC_discharge1'
                
                nextCtrlType = 'CV_discharge2';

              case 'CV_discharge2'
                
                nextCtrlType = 'CV_discharge2';

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
