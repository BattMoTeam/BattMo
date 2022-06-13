function nodeDependencyList = getNodeDependencyList(nodeList, A)
    nodeDependencyList = [];
    for inode = 1 : numel(nodeList)
        node = nodeList(inode);
        nodeDependencyList = vertcat(nodeDependencyList, node);
        parentnodes = find(A(:, node)');
        if any(parentnodes)
            nodeDependencyList = vertcat(nodeDependencyList, getNodeDependencyList(parentnodes, A));
        end
    end        
end


