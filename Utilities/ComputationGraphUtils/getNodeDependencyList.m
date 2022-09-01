function nodeDependencyList = getNodeDependencyList(nodeList, A, varargin)

    opt = struct('oneParentOnly', false);
    opt = merge_options(opt, varargin{:});
    
    nodeDependencyList = [];
    for inode = 1 : numel(nodeList)
        node = nodeList(inode);
        nodeDependencyList = vertcat(nodeDependencyList, node);
        parentnodes = find(A(:, node)');
        if any(parentnodes)
            if ~(opt.oneParentOnly)
                nodeDependencyList = vertcat(nodeDependencyList, getNodeDependencyList(parentnodes, A));
            else
                nodeDependencyList = vertcat(nodeDependencyList, parentnodes');
            end
        end
    end        
end


