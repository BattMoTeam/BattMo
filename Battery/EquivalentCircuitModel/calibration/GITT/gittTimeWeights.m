function weights = gittTimeWeights(time, groups)
%GITTTIMEWEIGHTS Trapezoidal quadrature weights for non-uniform samples.
%
% With groups, every nonempty group receives equal total weight. This lets
% current-on and current-off phases contribute equally without allowing a
% densely sampled phase to dominate the calibration score.

    time = double(time(:));
    assert(numel(time) >= 2 && all(diff(time) > 0), ...
           'Time must contain at least two strictly increasing samples.');
    dt = diff(time);
    weights = zeros(size(time));
    weights(1) = dt(1)/2;
    weights(end) = dt(end)/2;
    if numel(time) > 2
        weights(2:end - 1) = (dt(1:end - 1) + dt(2:end))/2;
    end

    if nargin < 2 || isempty(groups)
        weights = weights/sum(weights);
        return;
    end

    groups = string(groups(:));
    assert(numel(groups) == numel(time), 'groups and time must have equal lengths.');
    uniqueGroups = unique(groups, 'stable');
    totalPerGroup = 1/numel(uniqueGroups);
    for i = 1:numel(uniqueGroups)
        mask = groups == uniqueGroups(i);
        groupWeight = sum(weights(mask));
        assert(groupWeight > 0, 'A time-weight group has zero duration.');
        weights(mask) = totalPerGroup*weights(mask)/groupWeight;
    end
end
