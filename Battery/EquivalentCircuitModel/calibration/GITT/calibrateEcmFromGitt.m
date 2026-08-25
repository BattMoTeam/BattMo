function result = calibrateEcmFromGitt(dataSource, options)
%CALIBRATEECMFROMGITT Calibrate directional BattMo 2-RC models from GITT.
%
% SYNOPSIS:
%   result = calibrateEcmFromGitt(filename)
%   result = calibrateEcmFromGitt(table, options)
%
% The function has no file-writing side effects. It returns separate charge
% and discharge BattMo jsonstructs because the current ECM has no runtime
% mechanism for directional parameter switching. The 2-RC fits use
% unitBoxBFGS with gradients propagated by MRST ADI.

    if nargin < 2
        options = struct();
    end
    defaults = struct('nominalCapacity', 5, ...
                      'currentThreshold', 0.1, ...
                      'socTargets', [0.1, 0.3, 0.5, 0.7, 0.9], ...
                      'socGrid', (0:0.01:1)', ...
                      'runFullCalibration', false, ...
                      'maximumOcvIterations', 3, ...
                      'ocvConvergenceTolerance', 1e-4, ...
                      'preRestDuration', 60, ...
                      'postRestDuration', 1800, ...
                      'maximumFitSamples', 400, ...
                      'maximumOptimizerIterations', 60, ...
                      'multiStarts', 2, ...
                      'minimumTauRatio', 2, ...
                      'legacyReferenceFile', '');
    options = applyDefaults(options, defaults);
    options.socGrid = unique(double(options.socGrid(:)));
    assert(numel(options.socGrid) >= 2 && options.socGrid(1) == 0 && ...
           options.socGrid(end) == 1, 'socGrid must span SOC 0 through 1.');

    if istable(dataSource) && all(ismember({'Time_s', 'Voltage_V', 'Current_A', ...
                                            'Mode', 'Direction', 'Step', 'SOC'}, ...
                                           dataSource.Properties.VariableNames))
        data = dataSource;
    else
        data = loadGittData(dataSource, ...
                            'nominalCapacity', options.nominalCapacity, ...
                            'currentThreshold', options.currentThreshold);
    end
    segments = segmentGittData(data, ...
                               'preRestDuration', options.preRestDuration, ...
                               'postRestDuration', options.postRestDuration);
    ocvFits = estimateGittOcv(segments);
    ocvMaps = buildOcvMaps(ocvFits, options.socGrid);

    [fitIndices, validationIndices] = chooseSegments(segments, options);
    acceptedFits = struct([]);
    allFits = struct([]);
    previousFits = struct('Charge', [], 'Discharge', []);
    iterationMultiStarts = options.multiStarts;

    for iteration = 1:options.maximumOcvIterations
        iterationFits = struct([]);
        for directionCell = {'Charge', 'Discharge'}
            direction = directionCell{1};
            directionIndices = fitIndices(strcmp({segments(fitIndices).direction}, direction));
            [~, order] = sort([segments(directionIndices).soc]);
            directionIndices = directionIndices(order);
            previous = previousFits.(direction);
            for k = 1:numel(directionIndices)
                segment = segments(directionIndices(k));
                map = ocvMaps.(direction);
                fitOptions = {'maximumSamples', options.maximumFitSamples, ...
                              'minimumTauRatio', options.minimumTauRatio, ...
                              'maximumIterations', options.maximumOptimizerIterations, ...
                              'multiStarts', iterationMultiStarts, ...
                              'initialFit', previous};
                [fit, ~] = fitGittEcmSegment(segment, map.soc, map.voltage, fitOptions{:});
                fit.iteration = iteration;
                fit.segmentIndex = directionIndices(k);
                iterationFits = appendStruct(iterationFits, fit);
                if fit.accepted
                    previous = fit;
                end
            end
            previousFits.(direction) = previous;
        end
        allFits = [allFits, iterationFits]; %#ok<AGROW>
        acceptedFits = iterationFits([iterationFits.accepted]);
        if iteration >= options.maximumOcvIterations || isempty(acceptedFits)
            break;
        end
        [ocvMaps, maximumCorrection] = applyOcvCorrection(ocvMaps, acceptedFits, options.socGrid);
        if maximumCorrection < options.ocvConvergenceTolerance
            break;
        end
        iterationMultiStarts = 1;
    end

    for directionCell = {'Charge', 'Discharge'}
        direction = directionCell{1};
        directionFits = acceptedFits(strcmp({acceptedFits.direction}, direction));
        assert(numel(directionFits) >= 2, ...
               'Fewer than two accepted %s fits; cannot build a parameter map.', ...
               lower(direction));
        parameterMap = buildParameterMap(directionFits, options.socGrid);
        jsonstruct = buildBattMoJsonstruct(parameterMap, ocvMaps.(direction), ...
                                           options.nominalCapacity, direction);
        result.(lower(direction)).jsonstruct = jsonstruct;
        result.(lower(direction)).parameterMap = parameterMap;
        result.(lower(direction)).ocvMap = ocvMaps.(direction);
    end

    legacy = loadLegacyReference(options.legacyReferenceFile, dataSource);
    validation = validateSegments(segments, validationIndices, result, legacy);

    result.localFits = fitStructToTable(allFits);
    result.finalFits = fitStructToTable(acceptedFits);
    result.ocvFits = ocvFits;
    result.validation = validation;
    result.fitSegmentIndices = fitIndices(:);
    result.validationSegmentIndices = validationIndices(:);
    result.options = options;
    result.optimizationMethod = '2-RC: unitBoxBFGS with MRST ADI gradients';
    result.currentConvention = 'Positive current discharges the cell';
end

function options = applyDefaults(options, defaults)
    assert(isstruct(options) && isscalar(options), 'options must be a scalar struct.');
    names = fieldnames(defaults);
    for i = 1:numel(names)
        if ~isfield(options, names{i}) || isempty(options.(names{i}))
            options.(names{i}) = defaults.(names{i});
        end
    end
end

function maps = buildOcvMaps(ocvFits, socGrid)
    for directionCell = {'Charge', 'Discharge'}
        direction = directionCell{1};
        mask = strcmp(ocvFits.Direction, direction);
        assert(sum(mask) >= 2, 'Insufficient %s relaxation data for OCV.', lower(direction));
        [soc, value] = mergeDuplicateSoc(ocvFits.SOC(mask), ocvFits.OCV_V(mask));
        maps.(direction) = struct('soc', socGrid, ...
                                  'voltage', boundedPchip(soc, value, socGrid), ...
                                  'observedSocRange', [min(soc), max(soc)]);
    end
end

function [fitIndices, validationIndices] = chooseSegments(segments, options)
    allIndices = 1:numel(segments);
    if options.runFullCalibration
        fitIndices = allIndices;
        validationIndices = [];
        return;
    end
    fitIndices = [];
    validationIndices = [];
    for directionCell = {'Charge', 'Discharge'}
        direction = directionCell{1};
        candidates = find(strcmp({segments.direction}, direction));
        candidateSoc = [segments(candidates).soc];
        used = false(size(candidates));
        for target = options.socTargets(:)'
            distance = abs(candidateSoc - target);
            distance(used) = inf;
            [~, local] = min(distance);
            if isfinite(distance(local))
                fitIndices(end + 1) = candidates(local); %#ok<AGROW>
                used(local) = true;
                validationDistance = abs(candidateSoc - target);
                validationDistance(used) = inf;
                [~, validationLocal] = min(validationDistance);
                if isfinite(validationDistance(validationLocal))
                    validationIndices(end + 1) = candidates(validationLocal); %#ok<AGROW>
                    used(validationLocal) = true;
                end
            end
        end
    end
    fitIndices = unique(fitIndices, 'stable');
    validationIndices = unique(validationIndices, 'stable');
end

function [maps, maximumCorrection] = applyOcvCorrection(maps, fits, socGrid)
    maximumCorrection = 0;
    for directionCell = {'Charge', 'Discharge'}
        direction = directionCell{1};
        local = fits(strcmp({fits.direction}, direction));
        if numel(local) < 2
            continue;
        end
        [soc, correction] = mergeDuplicateSoc([local.soc]', [local.ocvOffset]');
        correction = max(-0.02, min(0.02, correction));
        denseCorrection = boundedPchip(soc, correction, socGrid);
        maps.(direction).voltage = maps.(direction).voltage + denseCorrection;
        maximumCorrection = max(maximumCorrection, max(abs(denseCorrection)));
    end
end

function map = buildParameterMap(fits, socGrid)
    [soc, R0] = mergeDuplicateSoc([fits.soc]', [fits.R0]');
    [~, R1] = mergeDuplicateSoc([fits.soc]', [fits.R1]');
    [~, tau1] = mergeDuplicateSoc([fits.soc]', [fits.tau1]');
    [~, R2] = mergeDuplicateSoc([fits.soc]', [fits.R2]');
    [~, tau2] = mergeDuplicateSoc([fits.soc]', [fits.tau2]');
    map.SOC = socGrid;
    map.R0 = exp(boundedPchip(soc, log(R0), socGrid));
    map.R1 = exp(boundedPchip(soc, log(R1), socGrid));
    map.tau1 = exp(boundedPchip(soc, log(tau1), socGrid));
    map.R2 = exp(boundedPchip(soc, log(R2), socGrid));
    map.tau2 = exp(boundedPchip(soc, log(tau2), socGrid));
    map.tau2 = max(map.tau2, 2*map.tau1);
    map.C1 = map.tau1./map.R1;
    map.C2 = map.tau2./map.R2;
    map.observedSocRange = [min(soc), max(soc)];
end

function jsonstruct = buildBattMoJsonstruct(map, ocv, nominalCapacity, direction)
    jsonstruct = struct();
    jsonstruct.nominalcellcapacity = nominalCapacity;
    if strcmp(direction, 'Charge')
        jsonstruct.initSOC = 0;
    else
        jsonstruct.initSOC = 1;
    end
    jsonstruct.lowerVoltageCutoff = min(ocv.voltage);
    jsonstruct.totalTime = 1;
    jsonstruct.I = tabulatedFunction([0, 1], [0, 0], {'time'});
    jsonstruct.OCP = tabulatedFunction(ocv.soc, ocv.voltage, {'SOC'});
    jsonstruct.R0 = tabulatedFunction(map.SOC, map.R0, {'SOC'});
    jsonstruct.R1 = tabulatedFunction(map.SOC, map.R1, {'SOC'});
    jsonstruct.C1 = tabulatedFunction(map.SOC, map.C1, {'SOC'});
    jsonstruct.R2 = tabulatedFunction(map.SOC, map.R2, {'SOC'});
    jsonstruct.C2 = tabulatedFunction(map.SOC, map.C2, {'SOC'});
end

function value = tabulatedFunction(x, y, arguments)
    value = struct('functionFormat', 'tabulated', ...
                   'argumentList', {arguments}, ...
                   'dataX', x(:)', ...
                   'dataY', y(:)');
end

function output = boundedPchip(x, y, query)
    [x, order] = sort(x(:));
    y = y(order);
    output = zeros(size(query));
    inside = query >= x(1) & query <= x(end);
    output(inside) = interp1(x, y, query(inside), 'pchip');
    output(query < x(1)) = y(1);
    output(query > x(end)) = y(end);
end

function [soc, value] = mergeDuplicateSoc(soc, value)
    [soc, ~, group] = unique(soc(:));
    value = accumarray(group, value(:), [], @median);
    [soc, order] = sort(soc);
    value = value(order);
end

function legacy = loadLegacyReference(filename, dataSource)
    legacy = [];
    if isempty(filename) && (ischar(dataSource) || isstring(dataSource))
        candidate = fullfile(fileparts(char(dataSource)), 'gitt_ecm_parameters.csv');
        if isfile(candidate)
            filename = candidate;
        end
    end
    if ~isempty(filename) && isfile(filename)
        legacy = readtable(filename, 'VariableNamingRule', 'preserve');
    end
end

function validation = validateSegments(segments, indices, result, legacy)
    n = numel(indices);
    pulseId = zeros(n, 1);
    soc = zeros(n, 1);
    direction = strings(n, 1);
    rmse2RC = zeros(n, 1);
    rmseLegacy1RC = nan(n, 1);
    for i = 1:n
        segment = segments(indices(i));
        directionName = segment.direction;
        configured = result.(lower(directionName));
        map = configured.parameterMap;
        ocvMap = configured.ocvMap;
        parameters = parametersAtSoc(map, segment.soc);
        parameters.U10 = 0;
        parameters.U20 = 0;
        parameters.ocvOffset = 0;
        parameters.ocvSlope = 0;
        ocv = interp1(ocvMap.soc, ocvMap.voltage, segment.sampleSoc, 'linear');
        predicted = simulateGittEcm(parameters, segment.time, segment.current, ...
                                    segment.sampleSoc, ocv);
        pulseId(i) = segment.pulseId;
        soc(i) = segment.soc;
        direction(i) = string(directionName);
        weights = gittTimeWeights(segment.time);
        rmse2RC(i) = sqrt(sum(weights.*(segment.voltage - predicted).^2));
        if ~isempty(legacy)
            rmseLegacy1RC(i) = legacyRmse(segment, legacy, directionName);
        end
    end
    validation = table(pulseId, soc, direction, rmse2RC, rmseLegacy1RC, ...
                       'VariableNames', {'PulseId', 'SOC', 'Direction', ...
                                         'RMSE2RC_V', 'RMSELegacy1RC_V'});
end

function parameters = parametersAtSoc(map, soc)
    parameters = struct('R0', interp1(map.SOC, map.R0, soc), ...
                        'R1', interp1(map.SOC, map.R1, soc), ...
                        'tau1', interp1(map.SOC, map.tau1, soc), ...
                        'R2', interp1(map.SOC, map.R2, soc), ...
                        'tau2', interp1(map.SOC, map.tau2, soc));
end

function rmse = legacyRmse(segment, legacy, direction)
    directionColumn = string(legacy.('Direction / 1'));
    subset = legacy(directionColumn == string(direction), :);
    x = subset.('State of Charge / 1');
    parameters = struct('R0', interp1(x, subset.('R0 / Ohm'), segment.soc), ...
                        'R1', interp1(x, subset.('R1 / Ohm'), segment.soc), ...
                        'tau1', interp1(x, subset.('Tau1 / s'), segment.soc), ...
                        'R2', 1e-8, 'tau2', 1);
    parameters.U10 = 0;
    parameters.U20 = 0;
    parameters.ocvOffset = 0;
    parameters.ocvSlope = 0;
    ocv = interp1(x, subset.('Open-circuit voltage / V'), segment.sampleSoc);
    predicted = simulateGittEcm(parameters, segment.time, segment.current, ...
                                segment.sampleSoc, ocv);
    weights = gittTimeWeights(segment.time);
    rmse = sqrt(sum(weights.*(segment.voltage - predicted).^2));
end

function tableValue = fitStructToTable(fits)
    if isempty(fits)
        tableValue = table();
        return;
    end
    pulseId = [fits.pulseId]';
    iteration = [fits.iteration]';
    direction = string({fits.direction})';
    soc = [fits.soc]';
    R0 = [fits.R0]';
    R1 = [fits.R1]';
    tau1 = [fits.tau1]';
    C1 = tau1./R1;
    R2 = [fits.R2]';
    tau2 = [fits.tau2]';
    C2 = tau2./R2;
    U10 = [fits.U10]';
    U20 = [fits.U20]';
    ocvOffset = [fits.ocvOffset]';
    ocvSlope = [fits.ocvSlope]';
    rmse = [fits.rmse]';
    maximumError = [fits.maximumError]';
    identifiabilityRatio = [fits.identifiabilityRatio]';
    boundLimited = [fits.boundLimited]';
    accepted = [fits.accepted]';
    tableValue = table(pulseId, iteration, direction, soc, R0, R1, tau1, C1, ...
                       R2, tau2, C2, U10, U20, ocvOffset, ocvSlope, rmse, ...
                       maximumError, identifiabilityRatio, boundLimited, accepted, ...
                       'VariableNames', {'PulseId', 'Iteration', 'Direction', 'SOC', ...
                       'R0_Ohm', 'R1_Ohm', 'Tau1_s', 'C1_F', 'R2_Ohm', 'Tau2_s', ...
                       'C2_F', 'U10_V', 'U20_V', 'OCVOffset_V', 'OCVSlope_V', ...
                       'RMSE_V', 'MaximumError_V', 'IdentifiabilityRatio', ...
                       'BoundLimited', 'Accepted'});
end

function array = appendStruct(array, value)
    if isempty(array)
        array = value;
    else
        array(end + 1) = value;
    end
end
