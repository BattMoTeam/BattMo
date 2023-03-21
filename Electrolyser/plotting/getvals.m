function [ymin, ymax, y] = getvals(str, states)
    str = sprintf('y = cellfun(@(state) state.%s, states, "uniformOutput", false);', str);
    eval(str);
    ymin = cellfun(@(v) min(v), y);
    ymax = cellfun(@(v) max(v), y);
end
    
