function data = loadGittData(dataSource, varargin)
%LOADGITTDATA Load GITT measurements into BattMo's canonical convention.
%
% The returned table contains Time_s, Voltage_V, Current_A, Temperature_C,
% Record, Mode, Direction, Step, Capacity_Ah and SOC. Current_A follows the
% BattMo convention: positive current discharges the cell.

    opt = struct('nominalCapacity', 5, ...
                 'currentThreshold', 0.1, ...
                 'temperatureTolerance', 3);
    opt = merge_options(opt, varargin{:});

    if istable(dataSource)
        raw = dataSource;
    elseif ischar(dataSource) || isstring(dataSource)
        filename = char(dataSource);
        assert(isfile(filename), 'GITT data file does not exist: %s', filename);
        [~, ~, extension] = fileparts(filename);
        switch lower(extension)
          case {'.csv', '.txt'}
            raw = readtable(filename, 'VariableNamingRule', 'preserve');
          otherwise
            error('Unsupported GITT data extension: %s. Convert measurements to CSV first.', ...
                  extension);
        end
    else
        error('dataSource must be a filename or a MATLAB table.');
    end

    names = raw.Properties.VariableNames;
    timeName = findColumn(names, {'Time_s', 'Test Time / s', 'Time_h', 'Test Time / h'});
    voltageName = findColumn(names, {'Voltage_V', 'Voltage / V', 'U_V'});
    currentName = findColumn(names, {'Current_A', 'Current / A', 'I_A'});
    temperatureName = findColumn(names, {'Temperature_C', 'Temperature / degC', 'T1_C'}, false);
    recordName = findColumn(names, {'Record', 'Record Count / 1', 'DataSet'}, false);

    time = double(raw.(timeName));
    if contains(lower(timeName), '_h') || contains(lower(timeName), '/ h')
        time = 3600*time;
    end
    voltage = double(raw.(voltageName));
    sourceCurrent = double(raw.(currentName));
    current = -sourceCurrent; % source: positive charge; BattMo: positive discharge

    if isempty(temperatureName)
        temperature = nan(size(time));
    else
        temperature = double(raw.(temperatureName));
    end
    if isempty(recordName)
        record = (1:numel(time))';
    else
        record = double(raw.(recordName));
    end

    valid = isfinite(time) & isfinite(voltage) & isfinite(current);
    time = time(valid);
    voltage = voltage(valid);
    current = current(valid);
    temperature = temperature(valid);
    record = record(valid);

    [time, order] = sort(time, 'ascend');
    voltage = voltage(order);
    current = current(order);
    temperature = temperature(order);
    record = record(order);
    [time, uniqueIndex] = unique(time, 'stable');
    voltage = voltage(uniqueIndex);
    current = current(uniqueIndex);
    temperature = temperature(uniqueIndex);
    record = record(uniqueIndex);

    assert(numel(time) >= 3, 'The GITT dataset has too few valid samples.');
    assert(all(diff(time) > 0), 'GITT time must be strictly increasing.');

    finiteTemperature = temperature(isfinite(temperature));
    if ~isempty(finiteTemperature)
        temperatureSpan = max(finiteTemperature) - min(finiteTemperature);
        assert(temperatureSpan <= opt.temperatureTolerance, ...
               'The dataset spans %.2f C; temperature surfaces are not supported.', ...
               temperatureSpan);
    end

    mode = repmat("Rest", numel(time), 1);
    mode(current < -opt.currentThreshold) = "Charge";
    mode(current > opt.currentThreshold) = "Discharge";
    step = ones(numel(time), 1);
    step(2:end) = 1 + cumsum(mode(2:end) ~= mode(1:end - 1));

    direction = mode;
    lastDirection = "";
    for i = 1:numel(direction)
        if mode(i) == "Rest"
            direction(i) = lastDirection;
        else
            lastDirection = mode(i);
        end
    end
    nextDirection = "";
    for i = numel(direction):-1:1
        if mode(i) == "Rest" && direction(i) == ""
            direction(i) = nextDirection;
        elseif mode(i) ~= "Rest"
            nextDirection = mode(i);
        end
    end

    chargeCurrent = -current;
    capacity = cumtrapz(time, chargeCurrent)/3600;
    capacitySpan = max(capacity) - min(capacity);
    assert(capacitySpan > 0, 'Cannot infer SOC because integrated capacity is constant.');
    soc = (capacity - min(capacity))/capacitySpan;

    relativeCapacityError = abs(capacitySpan - opt.nominalCapacity)/opt.nominalCapacity;
    if relativeCapacityError > 0.2
        warning('BattMo:GITT:CapacitySpan', ...
                'Integrated capacity span is %.3f Ah, versus nominal %.3f Ah.', ...
                capacitySpan, opt.nominalCapacity);
    end

    data = table(time, voltage, current, temperature, record, mode, direction, ...
                 step, capacity, soc, ...
                 'VariableNames', {'Time_s', 'Voltage_V', 'Current_A', ...
                                   'Temperature_C', 'Record', 'Mode', ...
                                   'Direction', 'Step', 'Capacity_Ah', 'SOC'});
end

function name = findColumn(names, candidates, required)
    if nargin < 3
        required = true;
    end
    name = '';
    for i = 1:numel(candidates)
        match = strcmpi(names, candidates{i});
        if any(match)
            name = names{find(match, 1)};
            return;
        end
    end
    if required
        error('Missing required data column. Expected one of: %s', ...
              strjoin(candidates, ', '));
    end
end
