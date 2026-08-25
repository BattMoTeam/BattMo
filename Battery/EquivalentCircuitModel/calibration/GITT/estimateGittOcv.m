function ocv = estimateGittOcv(segments, varargin)
%ESTIMATEGITTOCV Estimate directional OCV from GITT current-off relaxation.
%
% Each rest is represented as OCV + A1*exp(-t/tau1) +
% A2*exp(-t/tau2). Time constants are selected by deterministic grid
% refinement, while amplitudes and OCV are solved by weighted linear least
% squares. This initialization step does not require derivatives.

    opt = struct('maximumSamples', 160, ...
                 'coarseGridSize', 10, ...
                 'refinementLevels', 2, ...
                 'minimumTauRatio', 2, ...
                 'maximumTau', 36000, ...
                 'ocvMargin', 0.15);
    opt = merge_options(opt, varargin{:});

    n = numel(segments);
    pulseId = zeros(n, 1);
    soc = zeros(n, 1);
    direction = strings(n, 1);
    value = zeros(n, 1);
    stdev = zeros(n, 1);
    rmse = zeros(n, 1);
    tau1 = zeros(n, 1);
    tau2 = zeros(n, 1);
    usedFallback = false(n, 1);

    for i = 1:n
        segment = segments(i);
        indices = segment.postIndices;
        t = segment.time(indices);
        t = t - t(1);
        voltage = segment.voltage(indices);
        selected = evenlySpacedIndices(numel(t), opt.maximumSamples);
        t = t(selected);
        voltage = voltage(selected);

        [fit, ok] = fitDoubleExponential(t, voltage, opt);
        tailCount = max(3, ceil(0.1*numel(voltage)));
        tailValue = median(voltage(end - tailCount + 1:end));
        if ~ok || fit.ocv < min(voltage) - opt.ocvMargin || ...
                fit.ocv > max(voltage) + opt.ocvMargin
            fit.ocv = tailValue;
            fit.rmse = sqrt(sum(gittTimeWeights(t).*(voltage - tailValue).^2));
            fit.tau1 = nan;
            fit.tau2 = nan;
            usedFallback(i) = true;
        end

        pulseId(i) = segment.pulseId;
        soc(i) = mean(segment.sampleSoc(indices));
        direction(i) = string(segment.direction);
        value(i) = fit.ocv;
        rmse(i) = fit.rmse;
        stdev(i) = fit.rmse/sqrt(max(1, numel(voltage) - 3));
        tau1(i) = fit.tau1;
        tau2(i) = fit.tau2;
    end

    ocv = table(pulseId, soc, direction, value, stdev, rmse, tau1, tau2, ...
                usedFallback, ...
                'VariableNames', {'PulseId', 'SOC', 'Direction', 'OCV_V', ...
                                  'OCVStdev_V', 'RMSE_V', 'Tau1_s', ...
                                  'Tau2_s', 'UsedFallback'});
end

function [fit, ok] = fitDoubleExponential(time, voltage, opt)
    fit = struct('ocv', nan, 'rmse', inf, 'tau1', nan, 'tau2', nan);
    ok = false;
    if numel(time) < 5 || time(end) <= 0
        return;
    end

    positiveDt = diff(time);
    positiveDt = positiveDt(positiveDt > 0);
    if isempty(positiveDt)
        return;
    end
    tauMinimum = max(0.05, 2*median(positiveDt));
    tauMaximum = min(opt.maximumTau, max(4*tauMinimum, 2*time(end)));
    if tauMaximum <= opt.minimumTauRatio*tauMinimum
        return;
    end

    grid = logspace(log10(tauMinimum), log10(tauMaximum), opt.coarseGridSize);
    bestScore = inf;
    bestTau = [];
    bestCoefficients = [];
    for i = 1:numel(grid)
        for j = i + 1:numel(grid)
            if grid(j) < opt.minimumTauRatio*grid(i)
                continue;
            end
            [score, coefficients] = exponentialScore(grid(i), grid(j), time, voltage);
            if score < bestScore
                bestScore = score;
                bestTau = [grid(i), grid(j)];
                bestCoefficients = coefficients;
            end
        end
    end
    if isempty(bestTau)
        return;
    end

    logLower = log(tauMinimum);
    logUpper = log(tauMaximum);
    searchHalfWidth = (logUpper - logLower)/(opt.coarseGridSize - 1);
    for level = 1:opt.refinementLevels
        tau1Log = linspace(max(logLower, log(bestTau(1)) - searchHalfWidth), ...
                           min(logUpper - log(opt.minimumTauRatio), ...
                               log(bestTau(1)) + searchHalfWidth), ...
                           opt.coarseGridSize);
        tau2Log = linspace(max(logLower + log(opt.minimumTauRatio), ...
                               log(bestTau(2)) - searchHalfWidth), ...
                           min(logUpper, log(bestTau(2)) + searchHalfWidth), ...
                           opt.coarseGridSize);
        for i = 1:numel(tau1Log)
            for j = 1:numel(tau2Log)
                candidateTau = exp([tau1Log(i), tau2Log(j)]);
                if candidateTau(2) < opt.minimumTauRatio*candidateTau(1)
                    continue;
                end
                [candidateScore, candidateCoefficients] = exponentialScore( ...
                    candidateTau(1), candidateTau(2), time, voltage);
                if candidateScore < bestScore
                    bestScore = candidateScore;
                    bestTau = candidateTau;
                    bestCoefficients = candidateCoefficients;
                end
            end
        end
        searchHalfWidth = 2*searchHalfWidth/(opt.coarseGridSize - 1);
    end

    fit.ocv = bestCoefficients(1);
    fit.rmse = sqrt(bestScore);
    fit.tau1 = bestTau(1);
    fit.tau2 = bestTau(2);
    ok = isfinite(fit.ocv) && isfinite(fit.rmse);
end

function [score, coefficients] = exponentialScore(tau1, tau2, time, voltage)
    design = [ones(numel(time), 1), exp(-time/tau1), exp(-time/tau2)];
    weights = gittTimeWeights(time);
    rootWeights = sqrt(weights);
    coefficients = (rootWeights.*design)\(rootWeights.*voltage);
    residual = voltage - design*coefficients;
    score = sum(weights.*residual.^2);
end

function indices = evenlySpacedIndices(count, maximumCount)
    if count <= maximumCount
        indices = (1:count)';
    else
        indices = unique(round(linspace(1, count, maximumCount)))';
    end
end
