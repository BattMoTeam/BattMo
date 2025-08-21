classdef TimeControlModelInputParams < ControlModelInputParams

    properties

        usetable    % true if table is used
        usefunction % true if we use a matlab function

        %% property used in case usetable is true
        
        times        % Array with time value (should include at the end the end time, so that length(times) = length(durations) + 1)
        durations    % Array with duration value
        values       % Array with control value
        controltypes % Array with control type. The convention is
                     % - 1 for current
                     % - 2 for voltage

        
        %% property used in case usefunction is true

        controlValueFunction % function of time for the control value
        controlTypeFunction  % function of time for the type control value (1 : current, 2 : voltage)
        
    end

    methods
        
        function inputparams = TimeControlModelInputParams(jsonstruct);

            jsonstruct = setDefaultJsonStructField(jsonstruct, 'usetable', true);
            jsonstruct = setDefaultJsonStructField(jsonstruct, 'usefunction', false);

            usetable    = getJsonStructField(jsonstruct, 'usetable');
            usefunction = getJsonStructField(jsonstruct, 'usefunction');

            assert(xor(usetable, usefunction), 'options usetable and usefunction are not compatible');

            if usetable
                if isAssigned(jsonstruct, 'durations')
                    durations = getJsonStructField(jsonstruct, 'durations');
                    if isAssigned(jsonstruct, 'times')
                        error('cannot assign both durations and times')
                    else
                        jsonstruct.times = [0; cumsum(durations)];
                    end
                else
                    times = getJsonStructField(jsontruct, 'times');
                    jsonstruct.durations = diff(times);
                end
            end

            inputparams =  inputparams@ControlModelInputParams(jsonstruct);
            
        end
    end 
    
end

