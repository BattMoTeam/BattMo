classdef GenericControlModel < ControlModel

    properties
        
        controlsteps % cell array, each cell describing a control (see GenericControlModel.schema.json)
        
    end
    
    methods
        
        function model = GenericControlModel(inputparams)
            
            model = model@ControlModel(inputparams);
            
            fdnames = {'controlsteps'};
            model = dispatchParams(model, inputparams, fdnames);

            model = model.setupTerminationFunctions();
            
        end
        
        function model = registerVarAndPropfuncNames(model)
            
            model = registerVarAndPropfuncNames@ControlModel(model);
            
            varnames = {};
            % Control type : string that can take following value
            % - 'rest'
            % - 'current'
            % - 'voltage'
            varnames{end + 1} = 'ctrlType';
            % Control step index : integer that indicates the current control step in the controlsteps array
            varnames{end + 1} = 'ctrlStepIndex';
            % Time spanned by the current Control
            varnames{end + 1} = 'ctrlTime';
            % Time of the current control switch
            varnames{end + 1} = 'ctrlSwitchTime';
            % Time spanned since beginning
            varnames{end + 1} = 'time';
            
            model = model.registerVarNames(varnames);
            
            model = model.setAsStaticVarNames({'ctrlStepIndex', 'ctrlSwitchTime'});
            
            % Register the functions
            fn = @GenericControlModel.updateControlType;
            model = model.registerPropFunction({'ctrlType', fn, {'ctrlStepIndex'}});
            
            fn = @GenericControlModel.updateControlTime;
            model = model.registerPropFunction({'ctrlTime', fn, {'time', 'ctrlSwitchTime'}});

            fn = @GenericControlModel.updateControlEquation;
            model = model.registerPropFunction({'controlEquation', fn, {'ctrlType', 'E', 'I', 'ctrlTime'}});
            
        end
        
        function cleanState = addStaticVariables(model, cleanState, state)
            
            cleanState.ctrlStepIndex  = state.ctrlStepIndex;
            cleanState.ctrlSwitchTime = state.ctrlSwitchTime;

        end
        
        function newstate = addVariablesAfterConvergence(model, newstate, state)

            newstate.ctrlTime = state.ctrlTime;
            
        end
        
        function model = setupTerminationFunctions(model)

            controlsteps = model.controlsteps;

            for ictrl = 1 : numel(controlsteps)

                controlsteps{ictrl}.termination = GenericControlModel.setupComparisonFunction(controlsteps{ictrl}.termination);
            end

            model.controlsteps = controlsteps;
            
        end
        
        function schedule = setupSchedule(model, jsonstruct)

            schedule.is_used = false;

        end
        
        function [dt, done, currControl] = getTimeStep(model, itstep, schedule, state)
        % Returns the current time-step and control

            ctrlsteps = model.controlsteps;

            ctrlind   = state.ctrlStepIndex;
            
            if ctrlind > length(ctrlsteps)
                done = true;
                dt = [];
                currControl = [];
                return
            end

            dt = ctrlsteps{ctrlind}.timeStepSize;
            done = false;
            currControl = 1;
            
        end
        
        function state = updateControlType(model, state)
            
            istep    = model.getControlStep(state);
            ctrlstep = model.controlsteps{istep};
            
            state.ctrlType = ctrlstep.controltype;
            
        end

        function state = updateControlTime(model, state)

            state.ctrlTime = state.time - state.ctrlSwitchTime;
            
        end
        
        function [isCtrlDone, isLastControl] = isControlTerminated(model, state, state0)

            if ~isfield(state, 'ctrlType')
                state = model.updateControlType(state);
            end
            ctrlType  = state.ctrlType;
            [ictrlstep, isLastControl] = model.getControlStep(state);
            
            ctrlstep = model.controlsteps{ictrlstep};
            termination = ctrlstep.termination;
            
            isCtrlDone = termination.function(state);

        end
        
        function state = updateControlState(model, state, state0, dt)
            
        % This function is called by the solver at the end of a Newton step.  We check if we have triggered the
        % "termination" criteria from the control If the termination criteria is met, we increment the control step by
        % one, and proceed with the Newton algorithm
            
            [isCtrlDone, isLastControl] = model.isControlTerminated(state, state0);
            
            if isCtrlDone && ~isLastControl
                state.ctrlStepIndex  = state.ctrlStepIndex + 1;
                state.ctrlSwitchTime = 0;
                state = model.updateValueFromControl(state);
            end
            
        end

        
        function  [arefulfilled, state] = checkConstraints(model, state, state0, dt)
            
        % This function is called when the Newton method has converged, but we want to check if the solution we have
        % obtained does not break the termination criteria.
            
        end
        
        function doend = triggerSimulationEnd(model, state, state0_inner, drivingForces)

            [isCtrlDone, isLastControl] = model.isControlTerminated(state, state0_inner);
            
            doend = (isCtrlDone && isLastControl);
            
        end
        
        function state = updateValueFromControl(model, state)
        % From the given control type, set the corresponding control variable to the expected value
            
            ictrlstep   = model.getControlStep(state);
            controlstep = model.controlsteps{ictrlstep};
            ctrlType    = controlstep.controltype;
            
            state.ctrlType = ctrlType;
            
            switch ctrlType
                
              case 'rest'
                
                state.I = 0
                
              case 'current'
                
                givenI = controlstep.value;
                state.I = givenI;
                
              case 'voltage'
                
                givenE = controlstep.value;
                state.E = givenE;
                
              otherwise
                
                error('control type not recognized');
                
            end
            
        end

        function [ictrlstep, isLastControl] = getControlStep(model, state)
            
            ctrlsteps = model.controlsteps;
            ictrlstep = state.ctrlStepIndex;

            isLastControl = false;
            if ictrlstep == length(ctrlsteps)
                isLastControl = true;
            end

        end
        
        function state = updateControlEquation(model, state)
            
            ctrlType  = state.ctrlType;
            ictrlstep = state.ctrlStepIndex;
            
            controlstep = model.controlsteps{ictrlstep};
            
            switch ctrlType
                
              case 'rest'
                
                ctrleq = state.I;
                
              case 'current'
                
                givenI = controlstep.value;
                if strcmp(controlstep.direction, 'charge')
                    givenI = -givenI;
                end
                ctrleq = state.I - givenI;
                
              case 'voltage'
                
                givenE = controlstep.value;
                
                ctrleq = (state.E - givenE)*1e5;
                
              otherwise
                
                error('control type not recognized');
                
            end
            
            state.controlEquation = ctrleq;
            
        end
        
    end

    methods (Static)

        function isdone = basicComparison(val, refval, comparison)

            switch comparison

              case 'below'
                isdone = val < refval;
              case 'above'
                isdone = val > refval;
              case 'absolute value below'
                isdone = abs(val) < refval;
              case 'absolute value above'
                isdone = abs(val) > refval;
              otherwise
                error('comparison string not recognized')
                
            end

        end
        
        function termination = setupComparisonFunction(termination)
        % return a function that takes state as argument and returns true if the termination condition is fullfilled.

            comparison = termination.comparison;
            
            switch termination.quantity
              case 'current'
                tI = termination.value;
                termination.function = @(state) GenericControlModel.basicComparison(state.I, tI, comparison);
              case 'voltage'
                tE = termination.value;
                termination.function = @(state) GenericControlModel.basicComparison(state.E, tE, comparison);
              case 'time'
                tTime = termination.value;
                termination.function = @(state) GenericControlModel.basicComparison(state.ctrlTime, tTime, comparison);
              otherwise
                error('termination type not recognized')
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
