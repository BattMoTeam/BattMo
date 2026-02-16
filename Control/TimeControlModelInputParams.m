classdef TimeControlModelInputParams < ControlModelInputParams

    properties

        usetable % true if table is used

        times        % Array with time value (should include at the end the end time, so that length(times) = length(durations) + 1)
        durations    % Array with duration value
        values       % Array with control value
        controltypes % Array or scalar with control type. The convention is
                     % 1 for current
                     % 2 for voltage
                     % 3 for power

        usefunction % true if we use a matlab function

        functionname % function name, should be in matlab path

    end

    methods

        function inputparams = TimeControlModelInputParams(jsonstruct)

            jsonstruct = setDefaultJsonStructField(jsonstruct, 'usetable', true);
            jsonstruct = setDefaultJsonStructField(jsonstruct, 'usefunction', false);

            usetable    = getJsonStructField(jsonstruct, 'usetable');
            usefunction = getJsonStructField(jsonstruct, 'usefunction');

            assert(xor(usetable, usefunction), 'options usetable and usefunction are not compatible');

            if usetable
                if isAssigned(jsonstruct, 'filename')
                    filename = getJsonStructField(jsonstruct, 'filename');
                    s = load(filename);
                    jsonstruct.durations    = s.durations;
                    jsonstruct.times        = [0; cumsum(s.durations)];
                    jsonstruct.values       = s.values;
                    jsonstruct.controltypes = s.controltypes;
                else
                    if isAssigned(jsonstruct, 'durations')
                        durations = getJsonStructField(jsonstruct, 'durations');
                        if isAssigned(jsonstruct, 'times')
                            error('cannot assign both durations and times')
                        else
                            jsonstruct.times = [0; cumsum(durations)];
                        end
                    else
                        times = getJsonStructField(jsonstruct, 'times');
                        jsonstruct.durations = diff(times);
                    end
                end
            end

            controltypes = getJsonStructField(jsonstruct, 'controltypes');

            if isscalar(controltypes)
                jsonstruct.controltypes = controltypes * ones(size(jsonstruct.values));
            end

            inputparams =  inputparams@ControlModelInputParams(jsonstruct);

        end
    end

end
