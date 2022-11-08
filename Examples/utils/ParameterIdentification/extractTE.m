function [t, E, I] = extractTE(s)

    d1 = cellfun(@(x) [x.time; x.Control.E; x.Control.I], s, 'uniformoutput', false);
    d2 = [d1{:}]';
    t = d2(:, 1);
    E = d2(:, 2);
    I = d2(:, 3);

end