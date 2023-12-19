classdef CcCvControlModel < ControlModel
%
% Constant-Current-Constant-Voltage Control
%

    
    properties
        
        Imax

        lowerCutoffVoltage
        upperCutoffVoltage

        dIdtLimit
        dEdtLimit

        numberOfCycles

        % Control used initially. String that can take one of the following values
        % - 'discharging'
        % - 'charging'
        initialControl
    end
    
    
    methods

        function model = CcCvControlModel(inputparams)

            model = model@ControlModel(inputparams);
            
            fdnames = {'lowerCutoffVoltage', ...
                       'upperCutoffVoltage', ...
                       'dEdtLimit'         , ...
                       'dIdtLimit'         , ...
                       'numberOfCycles'    , ...
                       'initialControl'};
            
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
            % Variable to store number of cycles
            varnames{end + 1} = 'nCycles';
            model = model.registerVarNames(varnames);

            model = model.setAsStaticVarName('nCycles');
            
            fn = @CcCvControlModel.updateControlEquation;
            model = model.registerPropFunction({'controlEquation', fn, {'ctrlType', 'E', 'I'}});
            
        end

        function state = prepareStepControl(model, state, state0, dt, drivingForces)
            state.ctrlType = state0.nextCtrlType;
        end
        
        function state = updateControlEquation(model, state)
            
            Imax = model.Imax;
            Emin = model.lowerCutoffVoltage;
            Emax = model.upperCutoffVoltage;
            
            E = state.E;
            I = state.I;            
            ctrlType = state.ctrlType;
            
            switch ctrlType
              case 'CC_discharge1'
                ctrleq = I - Imax;
              case 'CC_discharge2'
                ctrleq = I;
              case 'CC_charge1'
                ctrleq = I + Imax;
              case 'CV_charge2'
                % TODO fix scaling
                ctrleq = (E - Emax)*1e5;
              otherwise
                error('ctrlType not recognized');
            end
            
            state.controlEquation = ctrleq;
            
        end

        function state = updateControlState(model, state)
            
            Emin = model.lowerCutoffVoltage;
            Emax = model.upperCutoffVoltage;
            Imax = model.Imax;
            
            ctrlType = state.ctrlType;
            E = state.E;
            I = state.I;
            
            if strcmp(ctrlType, 'CC_discharge1')
                if E <= Emin
                    state.ctrlType = 'CC_discharge2';
                    fprintf('switch control from CC_discharge1 to CC_discharge2\n');
                end
            end
            
            if strcmp(ctrlType, 'CV_discharge2')
                if I > Imax
                    state.ctrlType = 'CC_discharge1';
                    fprintf('switch control from CV to CC\n');
                end
            end    
            
        end
        

        function cleanState = addStaticVariables(model, cleanState, state)

            cleanState.numberOfCycles = state.numberOfCycles;
            cleanState.ctrlType       = state.ctrlType;
            
        end
        
        function  [arefulfilled, state] = checkConstraints(model, state)

            Imax    = model.Imax;
            Emin    = model.lowerCutoffVoltage;
            Emax    = model.upperCutoffVoltage;
            
            E        = state.E;
            I        = state.I;
            ctrlType = state.ctrlType;

            arefulfilled = true;
            switch ctrlType
              case 'CC_discharge1'
                if E < Emin
                    arefulfilled = false;
                    state.ctrlType = 'CC_discharge2';
                end
              case 'CC_discharge2'
                % do not check anything in this case
              case 'CC_charge1'
                if E > Emax
                    arefulfilled = false;
                    state.ctrlType = 'CV_charge2';
                end
              case 'CV_charge2'
                % do not check anything in this case
              otherwise
                
                error('controlType not recognized');
                
            end            
        end
        
        function state = updateControlAfterConvergence(model, state, state0, dt)

            Imax     = model.Imax;
            Emin     = model.lowerCutoffVoltage;
            Emax     = model.upperCutoffVoltage;
            dEdtMin  = model.dEdtLimit;
            dIdtMin  = model.dIdtLimit;
            initctrl = model.initialControl;
            E        = state.E;
            ctrlType = state.ctrlType;
            ncycles  = state0.numberOfCycles;
            
            dEdt = (state.E - state0.E)/dt;
            dIdt = (state.I - state0.I)/dt;
            
            switch ctrlType
              
              case 'CC_discharge1'
                
                nextCtrlType = 'CC_discharge1';
                if (E <= Emin) 
                    nextCtrlType = 'CC_discharge2';
                end
            
              case 'CC_discharge2'
                
                nextCtrlType = 'CC_discharge2';
                if (abs(dEdt) <= dEdtMin)
                    nextCtrlType = 'CC_charge1';
                    if strcmp(initctrl, 'charging')
                        ncycles = ncycles + 1;
                    end
                end
            
              case 'CC_charge1'

                nextCtrlType = 'CC_charge1';
                if (E >= Emax) 
                    nextCtrlType = 'CV_charge2';
                end 
                
              case 'CV_charge2'

                nextCtrlType = 'CV_charge2';
                if (abs(dIdt) < dIdtMin)
                    nextCtrlType = 'CC_discharge1';
                    if strcmp(initctrl, 'discharging')
                        ncycles = ncycles + 1;
                    end
                end                  
                
              otherwise
                
                error('controlType not recognized');
                
            end

            if ~strcmp(ctrlType, nextCtrlType)
                fprintf('Switch control type from %s to %s\n', ctrlType, nextCtrlType);
            end
            
            state.nextCtrlType   = nextCtrlType;
            state.numberOfCycles = ncycles;
            
        end

        function func = setupStopFunction(model)
            
            func = @(mainModel, state, state_prev) (state.Control.numberOfCycles >= mainModel.Control.numberOfCycles);
            
        end

        function func = setupControlFunction(model)
        % There is no control function in the sense that all the switches are automatically turned on and off.
            
            func = [];
            
        end

        function step = setupScheduleStep(model, timeSteppingParams)
            
        % Setup and a return the step structure that is part of the schedule which is used as input for
        % :mrst:`simulateScheduleAD`. For some control type, there is a natural construction for this structure. This is
        % why we include this method here, for convenience. It can be overloaded by derived classes. The
        % timeSteppingParams structure by default is given by the data described in :battmofile:`Utilities/JsonSchemas/TimeStepping.schema.json`

            if (nargin > 1)
                params = timeSteppingParams;
            else
                params = [];
            end

            params = model.parseTimeSteppingStruct(params);
            
            CRate   = model.CRate;
            ncycles = model.numberOfCycles;

            if ~isempty(params.totalTime) & isempty(ncycles)
                totalTime = params.totalTimes;
            else
                if ~isempty(ncycles)
                    warning('Both the total time and the number of cycles are given.\nWe do not use the given total time value but compute it instead from the number of cycles.');
                end
                totalTime = 2*ncycles*1.1*(1*hour/CRate);
            end

            if ~isempty(params.timeStepDuration)
                dt = params.timeStepDuration
            else
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
            control.CCCV = true;
            
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
