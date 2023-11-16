function inds = regexpSelect(names, name)
% implementation of regexp where whitespace matches every thing and name can be a cell in which case the match of each element will be included

    if iscell(name)
        inds = [];
        for iname = 1 : numel(name)
            addedinds = regexpSelect(names, name{iname});
            inds = vertcat(inds, addedinds);
        end
        inds = unique(inds);
    else
        name = regexprep(name, ' +', '.*');
        inds = regexp(names, name, 'once');
        inds = cellfun(@(x) ~isempty(x), inds);
        inds = find(inds);
    end
    
end
