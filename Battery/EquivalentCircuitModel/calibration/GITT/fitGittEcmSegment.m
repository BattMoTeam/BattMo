function [fit, prediction] = fitGittEcmSegment(segment, ocvSoc, ocvVoltage, varargin)
%FITGITTECMSEGMENT Fit a constrained 2-RC ECM to one GITT pulse window.

    opt = struct('maximumSamples', 400, ...
                 'voltageScale', 5e-3, ...
                 'huberThreshold', 2, ...
                 'minimumTauRatio', 2, ...
                 'maximumIterations', 60, ...
                 'multiStarts', 2, ...
                 'minimumIdentifiabilityRatio', 1e-10, ...
                 'initialFit', []);
    opt = merge_options(opt, varargin{:});

    selected = selectSamples(segment, opt.maximumSamples);
    time = segment.time(selected);
    voltage = segment.voltage(selected);
    current = segment.current(selected);
    sampleSoc = segment.sampleSoc(selected);
    mode = segment.mode(selected);
    ocv = interp1(ocvSoc, ocvVoltage, sampleSoc, 'pchip');

    positiveDt = diff(time);
    positiveDt = positiveDt(positiveDt > 0);
    assert(~isempty(positiveDt), 'The fitting window has no positive time steps.');
    tau1Minimum = max(0.05, 2*median(positiveDt));
    tau1Maximum = max(20, min(2000, segment.time(end)));
    tau2Minimum = max(1, opt.minimumTauRatio*tau1Minimum);
    tau2Maximum = max(100, min(36000, 2*segment.time(end)));
    assert(tau2Maximum > opt.minimumTauRatio*tau1Minimum, ...
           'The rest window cannot resolve two separated time constants.');

    lower = [-5, -5, log10(tau1Minimum), -5, log10(tau2Minimum), ...
             -0.3, -0.3, -0.15, -0.2]';
    upper = [log10(2), log10(2), log10(tau1Maximum), log10(2), ...
             log10(tau2Maximum), 0.3, 0.3, 0.15, 0.2]';

    initial = initialGuess(segment, lower, upper, opt.initialFit);
    starts = makeStarts(initial, opt.multiStarts);
    phase = repmat("Rest", numel(mode), 1);
    phase(mode ~= "Rest") = "Pulse";
    weights = gittTimeWeights(time, phase);
    objective = @(u) objectiveWithGradient(u, lower, upper, time, current, ...
                                           sampleSoc, ocv, voltage, weights, opt);

    delta1 = upper(3) - lower(3);
    delta2 = upper(5) - lower(5);
    inequality.A = zeros(1, 9);
    inequality.A(3) = delta1;
    inequality.A(5) = -delta2;
    inequality.b = -log10(opt.minimumTauRatio) - lower(3) + lower(5);

    bestValue = inf;
    bestU = initial;
    bestHistory = [];
    for i = 1:size(starts, 2)
        start = makeFeasible(starts(:, i), lower, upper, opt.minimumTauRatio);
        try
            [value, candidate, history] = unitBoxBFGS(start, objective, ...
                'maximize', false, ...
                'linIneq', inequality, ...
                'enforceFeasible', true, ...
                'maxIt', opt.maximumIterations, ...
                'objChangeTol', 1e-9, ...
                'gradTol', 1e-5, ...
                'lineSearchMaxIt', 20, ...
                'plotEvolution', false);
            if isfinite(value) && value < bestValue
                bestValue = value;
                bestU = candidate;
                bestHistory = history;
            end
        catch exception
            if isempty(exception.stack)
                location = 'unknown location';
            else
                location = sprintf('%s:%d', exception.stack(1).name, ...
                                   exception.stack(1).line);
            end
            warning('BattMo:GITT:OptimizationStart', ...
                    'Pulse %g optimization start %d failed at %s: %s', ...
                    segment.pulseId, i, location, exception.message);
        end
    end

    parameters = decodeParameters(bestU, lower, upper);
    fullOcv = interp1(ocvSoc, ocvVoltage, segment.sampleSoc, 'pchip');
    [prediction, ~, ~] = simulateGittEcm(parameters, segment.time, ...
                                         segment.current, segment.sampleSoc, fullOcv);
    residual = segment.voltage - prediction;
    rmse = sqrt(sum(gittTimeWeights(segment.time).*residual.^2));
    maximumError = max(abs(residual));

    jacobian = full(voltageJacobianAD(bestU, lower, upper, time, current, ...
                                      sampleSoc, ocv));
    rootWeights = sqrt(weights);
    weightedJacobian = rootWeights.*jacobian;
    singularValues = svd(weightedJacobian, 'econ');
    if isempty(singularValues) || singularValues(1) == 0
        identifiabilityRatio = 0;
    else
        identifiabilityRatio = singularValues(end)/singularValues(1);
    end
    mainAtBound = any(bestU(1:5) < 1e-6 | bestU(1:5) > 1 - 1e-6);

    degreesOfFreedom = max(1, numel(voltage) - numel(bestU));
    selectedResidual = voltage - simulateSelected(parameters, time, current, ...
                                                   sampleSoc, ocv);
    variance = sum(weights.*selectedResidual.^2)/degreesOfFreedom;
    covarianceU = variance*pinv(weightedJacobian'*weightedJacobian);
    physicalJacobian = physicalParameterJacobianAD(bestU, lower, upper);
    physicalCovariance = physicalJacobian*covarianceU*physicalJacobian';
    physicalStdev = sqrt(max(0, diag(physicalCovariance)));

    fit = parameters;
    fit.pulseId = segment.pulseId;
    fit.direction = segment.direction;
    fit.soc = segment.soc;
    fit.objective = bestValue;
    fit.rmse = rmse;
    fit.maximumError = maximumError;
    fit.identifiabilityRatio = identifiabilityRatio;
    fit.boundLimited = mainAtBound;
    fit.accepted = isfinite(bestValue) && all(isfinite(prediction)) && ...
                   parameters.tau2 >= opt.minimumTauRatio*parameters.tau1 && ...
                   ~mainAtBound && ...
                   identifiabilityRatio >= opt.minimumIdentifiabilityRatio;
    fit.stdevR0 = physicalStdev(1);
    fit.stdevR1 = physicalStdev(2);
    fit.stdevTau1 = physicalStdev(3);
    fit.stdevR2 = physicalStdev(4);
    fit.stdevTau2 = physicalStdev(5);
    fit.optimizerHistory = bestHistory;
    fit.normalizedParameters = bestU;
end

function [value, gradient] = objectiveWithGradient(u, lower, upper, time, current, ...
                                                   soc, ocv, measured, weights, opt)
    u = initVariablesADI(u);
    objective = objectiveValue(u, lower, upper, time, current, soc, ocv, ...
                               measured, weights, opt);
    value = objective.val;
    gradient = objective.jac{1}';
    if ~isfinite(value) || any(~isfinite(gradient))
        value = 1e12;
        gradient = zeros(size(gradient));
    end
end

function value = objectiveValue(u, lower, upper, time, current, soc, ocv, ...
                                measured, weights, opt)
    parameters = decodeParameters(u, lower, upper);
    predicted = simulateSelected(parameters, time, current, soc, ocv);
    scaledResidual = (measured - predicted)/opt.voltageScale;
    absoluteResidual = abs(scaledResidual);
    quadratic = absoluteResidual <= opt.huberThreshold;
    loss = 0*scaledResidual;
    if any(quadratic)
        loss(quadratic) = 0.5*scaledResidual(quadratic).^2;
    end
    if any(~quadratic)
        loss(~quadratic) = opt.huberThreshold*(absoluteResidual(~quadratic) - ...
                                               0.5*opt.huberThreshold);
    end
    value = sum(weights.*loss);
end

function voltage = simulateSelected(parameters, time, current, soc, ocv)
    voltage = simulateGittEcm(parameters, time, current, soc, ocv);
end

function parameters = decodeParameters(u, lower, upper)
    transformed = lower + (upper - lower).*u(:);
    parameters = struct('R0', 10.^transformed(1), ...
                        'R1', 10.^transformed(2), ...
                        'tau1', 10.^transformed(3), ...
                        'R2', 10.^transformed(4), ...
                        'tau2', 10.^transformed(5), ...
                        'U10', transformed(6), ...
                        'U20', transformed(7), ...
                        'ocvOffset', transformed(8), ...
                        'ocvSlope', transformed(9));
end

function initial = initialGuess(segment, lower, upper, previous)
    pulseStart = segment.pulseIndices(1);
    previousIndex = max(1, pulseStart - 1);
    deltaCurrent = segment.current(pulseStart) - segment.current(previousIndex);
    deltaVoltage = segment.voltage(pulseStart) - segment.voltage(previousIndex);
    if abs(deltaCurrent) > 1e-8
        R0 = min(0.5, max(1e-4, abs(deltaVoltage/deltaCurrent)));
    else
        R0 = 0.02;
    end
    pulseCurrent = median(abs(segment.current(segment.pulseIndices)));
    relaxationAmplitude = abs(segment.voltage(segment.postIndices(1)) - ...
                              segment.voltage(segment.postIndices(end)));
    dynamicResistance = relaxationAmplitude/max(pulseCurrent, 1e-3);
    R1 = max(1e-4, 0.6*dynamicResistance);
    R2 = max(1e-4, 0.4*dynamicResistance);
    physical = struct('R0', R0, 'R1', R1, 'tau1', 20, 'R2', R2, ...
                      'tau2', 200, 'U10', 0, 'U20', 0, ...
                      'ocvOffset', 0, 'ocvSlope', 0);
    if isstruct(previous) && ~isempty(previous)
        fields = fieldnames(physical);
        for i = 1:numel(fields)
            if isfield(previous, fields{i})
                physical.(fields{i}) = previous.(fields{i});
            end
        end
    end
    transformed = [log10(physical.R0), log10(physical.R1), log10(physical.tau1), ...
                   log10(physical.R2), log10(physical.tau2), physical.U10, ...
                   physical.U20, physical.ocvOffset, physical.ocvSlope]';
    initial = (transformed - lower)./(upper - lower);
    initial = min(0.98, max(0.02, initial));
end

function starts = makeStarts(initial, count)
    count = max(1, round(count));
    offsets = [zeros(9, 1), ...
               [0.08; -0.08; -0.12; 0.10; 0.12; 0; 0; 0; 0], ...
               [-0.08; 0.10; 0.10; -0.08; -0.10; 0; 0; 0; 0], ...
               [0; 0.12; -0.05; -0.12; 0.08; 0; 0; 0; 0]];
    starts = zeros(9, count);
    for i = 1:count
        starts(:, i) = min(0.98, max(0.02, initial + offsets(:, mod(i - 1, size(offsets, 2)) + 1)));
    end
end

function u = makeFeasible(u, lower, upper, minimumRatio)
    transformed = lower + (upper - lower).*u;
    if transformed(5) < transformed(3) + log10(minimumRatio)
        transformed(5) = min(upper(5), transformed(3) + log10(minimumRatio));
        if transformed(5) < transformed(3) + log10(minimumRatio)
            transformed(3) = transformed(5) - log10(minimumRatio);
        end
    end
    u = (transformed - lower)./(upper - lower);
    u = min(1, max(0, u));
end

function selected = selectSamples(segment, maximumSamples)
    n = numel(segment.time);
    if n <= maximumSamples
        selected = (1:n)';
        return;
    end
    transitions = find([true; segment.mode(2:end) ~= segment.mode(1:end - 1)]);
    preserve = [];
    for i = 1:numel(transitions)
        preserve = [preserve; (transitions(i) - 3:transitions(i) + 5)']; %#ok<AGROW>
    end
    preserve = preserve(preserve >= 1 & preserve <= n);
    remaining = max(0, maximumSamples - numel(unique(preserve)));
    if remaining > 0
        regular = round(linspace(1, n, remaining))';
    else
        regular = zeros(0, 1);
    end
    selected = unique([1; n; preserve; regular]);
end

function jacobian = voltageJacobianAD(u, lower, upper, time, current, soc, ocv)
    u = initVariablesADI(u);
    voltage = simulateSelected(decodeParameters(u, lower, upper), ...
                               time, current, soc, ocv);
    jacobian = voltage.jac{1};
end

function jacobian = physicalParameterJacobianAD(u, lower, upper)
    u = initVariablesADI(u);
    parameters = decodeParameters(u, lower, upper);
    physical = [parameters.R0; parameters.R1; parameters.tau1; ...
                parameters.R2; parameters.tau2];
    jacobian = physical.jac{1};
end
