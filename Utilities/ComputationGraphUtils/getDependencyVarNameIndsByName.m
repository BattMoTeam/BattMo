function depvarnaminds = getDependencyVarNameIndsByName(regstr, nodenames, A, varargin)
% return indices
    opt = struct('oneParentOnly', false);
    opt = merge_options(opt, varargin{:});
    
    varnameinds = regexp(nodenames, regstr, 'once');
    varnameinds = cellfun(@(x) ~isempty(x), varnameinds);
    varnameinds = find(varnameinds);
    depvarnaminds = getDependencyVarNameInds(varnameinds, A, 'oneParentOnly', opt.oneParentOnly);
    depvarnaminds = unique(depvarnaminds);
    
end

