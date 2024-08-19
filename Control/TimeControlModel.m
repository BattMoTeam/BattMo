classdef TimeControlModel < ControlModel


    properties

        usetable % true if table is used
        
        times        % Array with time value (should include at the end the end time, so that length(times) = length(durations) + 1)
        durations    % Array with time value
        values       % Array with control value
        controltypes % Array with control type. The convention is
                     % - 1 for current
                     % - 2 for voltage

        usefunction % true if we use a matlab function
        
        functionname % function name, should be in matlab path

        % Advanced parameters

        tolerance = 1e-4 % tolerance to skip timesteps (in second)
        
        % Helpers

        computeInput % function called to give update
         
        use_durations   % Setup when usetable is true
        functionhandler % Setup when usefunction is true
        
    end

    methods

        function model = TimeControlModel(inputparams)
            
            model = model@ControlModel(inputparams);

            fdnames = {'usetable'    , ...
                       'times'       , ...
                       'durations'   , ...
                       'values'      , ...
                       'controltypes', ...
                       'usefunction' , ...
                       'functionname' };
        
            model = dispatchParams(model, inputparams, fdnames);

            if model.usetable
                
                model.computeInput = @(t) model.computeInputFromTable(t);
                
            end

            if model.usefunction

                model.functionhandler = str2func(model.functionname);
                model.computeInput = @(t) model.computeInputFromFunction(t);
                
            end
            
        end
        
        function model = registerVarAndPropfuncNames(model)

            model = registerVarAndPropfuncNames@ControlModel(model);
            
            varnames = {};
            % Control type (string)
            % - 'constantCurrent'
            % - 'constantVoltage'
            varnames{end + 1} = 'ctrlType';            
            % control value that can be either a voltage or a current
            varnames{end + 1} = 'ctrlVal';            

            model = model.registerVarNames(varnames);
            
            fn = @CTimeControlModel.updateControlEquation;
            model = model.registerPropFunction({'controlEquation', fn, {'ctrlType', 'ctrlVal', 'E', 'I'}});
            
        end

        function [ctrlVal, ctrlType] = computeInputFromTable(model, t)

        % We could be more efficient here and keep track of previous index to avoid full search (the time spent for
        % that is probabely negligeable compared to the rest.)
            ind = find(t >= model.times, 1, 'last');

            if t > model.times(end)
                error('outside of time table')
            elseif t == model.times(end)
                ind = numel(model.values);
            end
            
            ctrlVal  = model.values(ind);
            ctrlType = model.controltypes(ind);

            switch ctrlType
              case 1
                ctrlType = 'constantCurrent';
              case 2
                ctrlType = 'constantVoltage';
              otherwise
                error('ctrlType not recognized. It should be equal to 1 or 2')
            end
            
        end
        
        function [ctrlVal, ctrlType] = computeInputFromFunction(model, t)

            [ctrlVal, ctrlType] = model.functionhandler(t);

        end

        function state = updateControlEquation(model, state)
            
            E        = state.E;
            I        = state.I;            
            ctrlVal  = state.ctrlVal;
            ctrlType = state.ctrlType;

            switch ctrlType
                
              case 'constantCurrent'
                
                ctrleq = I - ctrlVal;
                
              case 'constantVoltage'
                
                %% TODO : fix hard-coded scaling
                ctrleq = (E - ctrlVal)*1e5;
                
              otherwise
                
                error('ctrlType not recognized');
                
            end
            
            state.controlEquation = ctrleq;
            
        end
        
        function cleanState = addStaticVariables(model, cleanState, state)

            cleanState.ctrlType = state.ctrlType;
            
        end

        function func = setupControlFunction(model)

            func = [];
            
        end
        
        function step = setupScheduleStep(model, timeSteppingParams)
            
            if (nargin > 1)
                params = timeSteppingParams;
            else
                params = [];
            end

            if model.usetable

            end
            % Call parser for TimeStepping structure with some default values
            params = model.parseTimeSteppingStruct(params);

            totalTime = model.times(end);

            givendt = false;
            if ~isempty(params.timeStepDuration)
                dt = params.timeStepDuration;
                givendt = true;
            elseif ~isempty(params.timeStepDuration)
                n  = params.numberOfTimeSteps;
                dt = totalTime/n;
                givendt = true;
            end

            if givendt
                if params.useRampup
                    n = params.numberOfRampupSteps;
                else
                    n = 0;
                end
                
                dts = rampupTimesteps(totalTime, dt, n);

            end

            if model.usetable
                
                time = [0; cumsum(model.durations)];

                if givendt
                    time1 = [0; cumsum(dts)];
                    
                    time =  [time; time1];
                    time =  sort(time);
                    time = uniquetol(time, model.tolerance);
                end
                
                dts = diff(time);
                
            end
            
            step = struct('val', dts, 'control', ones(numel(dts), 1));

        end
        
    end
    
    
end

