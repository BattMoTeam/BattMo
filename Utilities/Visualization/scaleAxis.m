function scaleAxis(ax, dim, scaling)

    % dim is 'x', 'y', or 'z'

    dim = upper(dim);
    tick = get(ax, sprintf('%sTick', dim));
    tickvals = tick * scaling;
    ticklabel = arrayfun(@num2str, tickvals, 'un', false);
    set(ax, sprintf('%sTickLabel', dim), ticklabel);

end
