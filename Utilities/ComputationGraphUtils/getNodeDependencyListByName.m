function nodeDependencyList = getNodeDependencyListByName(regstr, nodenames, A)
    nodeList = regexp(nodenames, regstr, 'once');
    nodeList = cellfun(@(x) ~isempty(x), nodeList);
    nodeList = find(nodeList);
    nodeDependencyList = getNodeDependencyList(nodeList, A);
    nodeDependencyList = unique(nodeDependencyList);
end

