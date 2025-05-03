function jsonstruct = convertPyBaMMtoJson(experiment)

    if ischar(experiment)
        experiment = {experiment};
    end

    controlsteps = {};

    for iexp = 1:length(experiment)

        step = experiment{iexp};
        step_dict = struct();

        if containsi(step, 'Rest')

            step_dict.controltype = 'rest';
            [values, units] = extract_numeric_values(step);
            value = convert_to_seconds(values, units);

            step_dict.termination = struct('quantity', 'time', ...
                                           'value', value);
            step_dict.timeStepSize = 600;

        elseif containsi(step, 'Discharge') || containsi(step, 'Charge')

            if containsi(step, 'Discharge')
                direction  = 'discharge';
                comparison = 'below';
            else
                direction  = 'charge';
                comparison = 'above';
            end

            [values, units] = extract_numeric_values(step);
            [value1, quantity1] = type_to_unit(values(1), units{1});
            assert(strcmpi(quantity1, 'current'), 'Cannot %s with %s, can only use current', direction, quantity1);
            step_dict.controltype = quantity1;

            [value2, quantity2] = type_to_unit(values(2), units{2});

            step_dict.value = value1;
            step_dict.direction = direction;
            step_dict.termination = struct('quantity', quantity2, ...
                                           'value', value2, ...
                                           'comparison', comparison);
            step_dict.timeStepSize = 600;

        elseif containsi(step, 'Hold')

            [values, units] = extract_numeric_values(step);
            step_dict.controltype = 'voltage';
            value1 = convert_to_V(values(1), units{1});

            [value2, quantity2] = type_to_unit(values(2), units{2});
            step_dict.value = value1;
            step_dict.termination = struct('quantity', quantity2, ...
                                           'value', value2, ...
                                           'comparison', 'absolute value below');
            step_dict.timeStepSize = 600;

        else
            error('Unknown control step: %s', step);
        end

        controlsteps{end+1} = step_dict; %#ok<AGROW>
    end

    jsonstruct = struct('Control', ...
                        struct('controlPolicy', 'Generic', ...
                               'controlsteps', {controlsteps}));

end


function [values, units] = extract_numeric_values(str)
% Extract numeric values and their corresponding units
    tokens = regexp(str, '([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)\s*([a-zA-Z]*)', 'tokens');
    values = cellfun(@(x) str2double(x{1}), tokens);
    units = cellfun(@(x) x{2}, tokens, 'UniformOutput', false);

    if isscalar(values)
       units = units{1};
    end
end


function [value, quantity] = type_to_unit(value, unit)

    if containsi(unit, 'V')
        quantity = 'voltage';
        value = convert_to_V(value, unit);
    elseif containsi(unit, 'A')
        quantity = 'current';
        value = convert_to_A(value, unit);
    elseif containsi(unit, time_units())
        quantity = 'time';
        value = convert_to_seconds(value, unit);
    elseif containsi(unit, 'W')
        quantity = 'power';
        value = convert_to_W(value, unit);
    else
        error('Unknown unit: %s', unit);
    end

end


function [t, s, m, h] = time_units()

    s = struct('s', 1, ...
               'second', 1, ...
               'seconds', 1);

    m = struct('min', minute, ...
               'minute', minute, ...
               'minutes', minute, ...
               'mins', minute);

    h = struct('h', hour, ...
               'hour', hour, ...
               'hours', hour);

    s_fields = fieldnames(s);
    m_fields = fieldnames(m);
    h_fields = fieldnames(h);

    t = [s_fields; m_fields; h_fields];

end


function seconds = convert_to_seconds(value, unit)

    assert(isscalar(value), 'Only scalar values are supported');

    [t, s, m, h] = time_units();

    switch lower(unit)
      case fieldnames(s)
        seconds = value;
      case fieldnames(m)
        seconds = value * minute;
      case fieldnames(h)
        seconds = value * hour;
      otherwise
        error('Unknown unit: %s', unit);
    end

end


function V = convert_to_V(value, unit)

    assert(isscalar(value), 'Only scalar values are supported');

    switch unit
      case 'V'
        V = value;
      case 'mV'
        V = value*milli;
      otherwise
        error('Unknown unit: %s', unit);
    end

end


function A = convert_to_A(value, unit)

    assert(isscalar(value), 'Only scalar values are supported');

    switch unit
      case 'A'
        A = value;
      case 'mA'
        A = value*milli;
      otherwise
        error('Unknown unit: %s', unit);
    end

end

function W = convert_to_W(value, unit)

    assert(isscalar(value), 'Only scalar values are supported');

    switch unit
      case 'W'
        W = value;
      case 'mW'
        W = value*milli;
      otherwise
        error('Unknown unit: %s', unit);
    end

end


function s = containsi(a, b)

    s = contains(a, b, 'IgnoreCase', true);

end
