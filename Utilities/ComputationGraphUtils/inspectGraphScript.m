function inspectGraphScript(model)

    cgf = ComputationalGraphFilter(model);
    
    plotgraph = true;
    
    if plotgraph
        % cgf.includeNodeNames = 'cSurface';
        % cgf.includeNodeNames = 'phiElectrolyte';
        % gg = cgf.setupDescendantGraph();
        % gg = cgf.getFilteredGraph('oneParentOnly', true);    
        [g, edgelabels] = cgf.getFilteredGraph();
        figure
        h = plot(g, 'edgelabel', edgelabels, 'nodefontsize', 10);
    end


    A = cgf.A;
    nodenames = cgf.nodenames;

    doprintspecialvariables = true;

    if doprintspecialvariables
        fprintf('Root variables \n');
        nodenames(all(A == 0, 1))
        fprintf('Tail variables \n');
        nodenames(all(A' == 0, 1))
    end


    p = topological_order(A);


    funcCallList = {};
    for ind = 1 : numel(p)
        iprop = full(A(:, p(ind)));
        iprop = unique(iprop(iprop>0));
        if ~isempty(iprop)
            assert(numel(iprop) == 1, 'problem');
            propfunction = cgf.model.propertyFunctionList{iprop};
            fn = propfunction.fn;
            mn = propfunction.modelnamespace;
            mn = join(mn, '.');
            if ~isempty(mn)
                mn = mn{1};
                statename = sprintf('state.%s', mn);
            else
                statename = 'state';
            end
            fnname = func2str(fn);
            fnname = regexp(fnname, "\.(.*)", 'tokens');
            fnname = fnname{1}{1};
            fnname = horzcat(mn, {fnname});
            fnname = join(fnname, '.');
            fnname = fnname{1};

            funcCallList{end + 1} = sprintf('%s = model.%s(%s);', statename, fnname, statename);
        end
    end


    [~, ia, ic] = unique(funcCallList, 'first');

    ia = sort(ia);
    reducedFuncCallList = funcCallList(ia);

    fprintf('Function call list\n');
    for ind = 1 : numel(reducedFuncCallList)
        fprintf('%s\n', reducedFuncCallList{ind});
    end
    
end