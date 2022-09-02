function nodeDependencyList = getNodeDependencyListByName(regstr, nodenames, A, varargin)

    opt = struct('oneParentOnly', false);
    opt = merge_options(opt, varargin{:});
    
    nodeList = regexp(nodenames, regstr, 'once');
    nodeList = cellfun(@(x) ~isempty(x), nodeList);
    nodeList = find(nodeList);
    nodeDependencyList = getNodeDependencyList(nodeList, A, 'oneParentOnly', opt.oneParentOnly);
    nodeDependencyList = unique(nodeDependencyList);
    
end

