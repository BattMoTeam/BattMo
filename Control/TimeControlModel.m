classdef TimeControlModel < ControlModel


    properties

        inputtype % type used to give the time control input. It can be
                  % - table
                  % - function

        %% property used in case of table input
        
        times        % Array with time value (should include at the end the end time, so that length(times) = length(durations) + 1)
        durations    % Array with duration value
        values       % Array with control value
        controltypes % Array with control type. The convention is
                     % - 1 for current
                     % - 2 for voltage

        
        %% property used in case of function input

        controlValueFunction % function of time for the control value
        controlTypeFunction  % function of time for the type control value (1 : current, 2 : voltage)
        
        %% Helpers 

        usetable    % true if table is used
        usefunction % true if we use a matlab function

        computeInput % function called to give update
        
        use_durations   % Setup when usetable is true
        
        % Function handlers instantiated from controlValueFunction controlTypeFunction  

        controlValueFunc 
        controlTypeFunc
        
    end

    methods

        function model = TimeControlModel(inputparams)
            
            model = model@ControlModel(inputparams);

            fdnames = {'inputtype'           , ...
                       'times'               , ...
                       'durations'           , ...
                       'values'              , ...
                       'controltypes'        , ...
                       'controlValueFunction', ...
                       'controlTypeFunction'};

            model = dispatchParams(model, inputparams, fdnames);

            switch model.inputtype
              case 'table'
                model.usetable = true;
              case 'function'
                model.usefunction = true;
              otherwise
                error('input type not recognized');
            end
            
            if model.usetable
                
                model.computeInput = @(t) model.computeInputFromTable(t);
                
            end

            if model.usefunction

                model.controlValueFunc = setupFunction(model.controlValueFunction);
                model.controlTypeFunc  = setupFunction(model.controlTypeFunction);
                
                model.computeInput = @(t) model.computeInputFromFunction(t);

            end
            
        end
        
        function model = registerVarAndPropfuncNames(model)

            model = registerVarAndPropfuncNames@ControlModel(model);
            
            varnames = {};
            % Control type 
            % - 1 : Current control
            % - 2 : Voltage control
            varnames{end + 1} = 'ctrlType';            
            % control value that can be either a voltage or a current
            varnames{end + 1} = 'ctrlVal';            
            
            model = model.registerVarNames(varnames);
            
            fn = @CTimeControlModel.updateControlEquation;
            model = model.registerPropFunction({'controlEquation', fn, {'ctrlType', 'ctrlVal', 'E', 'I'}});
            
        end

        
        <<<<<<< HEAD
        function [ctrlVal, ctrlType] = computeInput(model, t)

            ctrlVal  = model.getValue.eval(t);
            ctrlType = model.getType.eval(t);
            =======
            function [ctrlVal, ctrlType] = computeInputFromFunction(model, t)

                ctrlVal  = model.controlValueFunc(t);
                ctrlType = model.controlTypeFunc(t);

                switch ctrlType
                  case 1
                    ctrlType = 'constantCurrent';
                  case 2
                    ctrlType = 'constantVoltage';
                  otherwise
                    error('ctrlType not recognized. It should be equal to 1 or 2')
                end
                
            end

            function state = updateControlEquation(model, state)
                
                E        = state.E;
                I        = state.I;            
                ctrlVal  = state.ctrlVal;
                ctrlType = state.ctrlType;

                switch ctrlType
                    
                  case 1
                    
                    ctrleq = I - ctrlVal;
                    
                  case 2
                    
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


                if model.usetable
                    totalTime = model.times(end);
                else
                    totalTime = params.totalTime;
                end

                if isa(totalTime, 'UnAssigned')
                    if isa(model.getValue, 'TabulatedFunction1D')
                        totalTime = model.getValue.dataX(end);
                    else
                        error('total time is not given and the input function is not a tabulated function');
                    end
                end
                
                if isAssigned(params, {'timeStepDuration'})
                    dt = params.timeStepDuration;
                    givendt = true;
                elseif ~isempty(params.numberOfTimeSteps)
                    n  = params.numberOfTimeSteps;
                    dt = totalTime/n;
                end

                dts = rampupTimesteps(totalTime, dt, 0);
                
                step = struct('val', dts, 'control', ones(numel(dts), 1));

            end
            
        end
        
        
    end



    
end

