function depvarnameinds = getDependencyVarNameInds(varnameinds, A, varargin)
% for given list of index varnameinds, returns a list with index of nodes (or v

    opt = struct('oneParentOnly', false);
    opt = merge_options(opt, varargin{:});
    
    depvarnameinds = [];
    for inode = 1 : numel(varnameinds)
        node = varnameinds(inode);
        depvarnameinds = vertcat(depvarnameinds, node);
        parentnodes = find(A(:, node)');
        if any(parentnodes)
            if ~(opt.oneParentOnly)
                depvarnameinds = vertcat(depvarnameinds, getDependencyVarNameInds(parentnodes, A));
            else
                depvarnameinds = vertcat(depvarnameinds, parentnodes');
            end
        end
    end        
end


