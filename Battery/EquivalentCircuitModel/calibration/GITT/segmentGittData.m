function segments = segmentGittData(data, varargin)
%SEGMENTGITTDATA Partition canonical GITT data into rest-pulse-rest windows.

    opt = struct('preRestDuration', 60, ...
                 'postRestDuration', 1800, ...
                 'minimumPulseSamples', 3, ...
                 'minimumRestSamples', 5);
    opt = merge_options(opt, varargin{:});

    required = {'Time_s', 'Voltage_V', 'Current_A', 'Mode', ...
                'Direction', 'Step', 'SOC'};
    assert(all(ismember(required, data.Properties.VariableNames)), ...
           'Input is not a canonical GITT table. Call loadGittData first.');

    steps = unique(data.Step, 'stable');
    segments = struct('pulseId', {}, 'direction', {}, 'soc', {}, ...
                      'time', {}, 'voltage', {}, 'current', {}, ...
                      'sampleSoc', {}, 'mode', {}, 'preIndices', {}, ...
                      'pulseIndices', {}, 'postIndices', {});

    for k = 2:numel(steps) - 1
        pulseStep = steps(k);
        pulseIndices = find(data.Step == pulseStep);
        if isempty(pulseIndices) || data.Mode(pulseIndices(1)) == "Rest"
            continue;
        end

        preIndices = find(data.Step == steps(k - 1));
        postIndices = find(data.Step == steps(k + 1));
        if isempty(preIndices) || isempty(postIndices) || ...
                data.Mode(preIndices(1)) ~= "Rest" || ...
                data.Mode(postIndices(1)) ~= "Rest"
            continue;
        end
        if numel(pulseIndices) < opt.minimumPulseSamples || ...
                numel(postIndices) < opt.minimumRestSamples
            continue;
        end

        pulseStart = data.Time_s(pulseIndices(1));
        pulseEnd = data.Time_s(pulseIndices(end));
        preIndices = preIndices(data.Time_s(preIndices) >= pulseStart - opt.preRestDuration);
        postIndices = postIndices(data.Time_s(postIndices) <= pulseEnd + opt.postRestDuration);
        if isempty(preIndices) || numel(postIndices) < opt.minimumRestSamples
            continue;
        end

        indices = [preIndices; pulseIndices; postIndices];
        segment = struct();
        segment.pulseId = double(pulseStep);
        segment.direction = char(data.Direction(pulseIndices(1)));
        segment.soc = mean(data.SOC(pulseIndices));
        segment.time = data.Time_s(indices) - data.Time_s(indices(1));
        segment.voltage = data.Voltage_V(indices);
        segment.current = data.Current_A(indices);
        segment.sampleSoc = data.SOC(indices);
        segment.mode = data.Mode(indices);
        segment.preIndices = (1:numel(preIndices))';
        segment.pulseIndices = (numel(preIndices) + (1:numel(pulseIndices)))';
        segment.postIndices = (numel(preIndices) + numel(pulseIndices) + ...
                               (1:numel(postIndices)))';
        segments(end + 1) = segment; %#ok<AGROW>
    end

    assert(~isempty(segments), 'No complete rest-pulse-rest GITT windows were found.');
end
