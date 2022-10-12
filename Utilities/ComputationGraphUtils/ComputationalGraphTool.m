classdef ComputationalGraphTool

    properties
        graph
        model
        A
        nodenames
        includeNodeNames
        excludeNodeNames
    end
    
    methods
        
        function cgt = ComputationalGraphTool(model)
            if isempty(model.propertyFunctionList)
                model = model.registerVarAndPropfuncNames();
            end
            cgt.model = model;
            g = setupGraph(model);
            cgt.graph = g;
            % In adjacency matrix A, 
            % - column index      : output variable index (as in model.varNameList)
            % - row index         : input variable (as in model.varNameList)
            % - coefficient value : property function index as in model.propertyFunctionList
            cgt.A         = adjacency(g, 'weighted');
            cgt.nodenames = g.Nodes.Variables;

            
        end

        function [propfuncs, nodenames] = findPropFunction(cgt, nodename)
            
            nodenames = cgt.nodenames;
            A         = cgt.A;
            model     = cgt.model;
            
            indSelectedNodenames = regexp(nodenames, nodename, 'once');
            indSelectedNodenames = cellfun(@(x) ~isempty(x), indSelectedNodenames);
            indPropfunctions = A(:, indSelectedNodenames);
            indPropfunctions = unique(indPropfunctions(:));

            propfuncs = model.propertyFunctionList(indPropfunctions(indPropfunctions > 0));

            if numel(propfuncs) == 1
                propfuncs = propfuncs{1};
            end
            
        end

        function printPropFunction(cgt, nodename)
            propfuncs = cgt.findPropFunction(nodename);

            if numel(propfuncs) == 1
                propfuncs = {propfuncs};
            end

            for iprop = 1 : numel(propfuncs)
                propfunc = propfuncs{iprop};
                varname = propfunc.varname;
                varname_s = varname.resolveIndex();
                for ind = 1 : numel(varname_s)
                    fullname = varname_s{ind}.getIndexedFieldname();
                    fprintf('state.%s <- ', fullname);
                end
                fprintf(' (%s) <- ', propfuncs{iprop}.getFunctionSignature());
                inputvarnames = propfunc.inputvarnames;
                inputvarstrs = {};
                for ind = 1 : numel(inputvarnames)
                    varname = inputvarnames{ind};
                    inputvarstrs{end + 1} = sprintf('state.%s', varname.getFieldname());
                end
                inputvarstr = join(inputvarstrs, {', '});
                fprintf('[%s]\n', inputvarstr{1});
            end
            
        end
        
        function [g, edgelabels] = getComputationalGraph(cgt, varargin)
            
            opt = struct('type', 'ascendant', ...
                         'oneParentOnly', false);
            opt = merge_options(opt, varargin{:});

            propfunctionlist = cgt.model.propertyFunctionList;
            nodenames        = cgt.nodenames;
            includeNodeNames = cgt.includeNodeNames;
            excludeNodeNames = cgt.excludeNodeNames;

            A = (cgt.A);
            if strcmp(opt.type, 'descendant')
                A = A';
            end
            
            if ~isempty(includeNodeNames)
                
                nodes = getNodeDependencyListByName(includeNodeNames, nodenames, A, 'oneParentOnly', opt.oneParentOnly);

                nodenames = nodenames(nodes);
                A = A(nodes, nodes);
                
                if ~isempty(excludeNodeNames)
                    removenodes = regexp(nodenames, excludeNodeNames, 'once');
                    removenodes = cellfun(@(x) ~isempty(x), removenodes);
                    nodenames = nodenames(~removenodes);
                    A = A(~removenodes, ~removenodes);
                end

            end
            
            if strcmp(opt.type, 'descendant')
                A = A';
            end
            
            g = digraph(A, nodenames);

            propinds = g.Edges.Weight;

            if nargout > 1
                edgelabels = {};
                for iprop = 1 : numel(propinds)
                    propind = propinds(iprop);
                    edgelabels{end + 1} = func2str(propfunctionlist{propind}.fn);
                end
            end
            
        end
        
        function [g, edgelabels] = setupAscendantGraph(cgt, varargin)
            [g, edgelabels] = cgt.setupGraph('type', 'ascendant', varargin{:});
        end
        
        function [g, edgelabels] = setupDescendantGraph(cgt, varargin)
            [g, edgelabels] = cgt.setupGraph('type', 'descendant', varargin{:});            
        end


        function printSpecialVariables(cgt)

            A = cgt.A;
            nodenames = cgt.nodenames;
            
            fprintf('Root variables \n');
            nodenames(all(A == 0, 1))
            fprintf('Tail variables \n');
            nodenames(all(A' == 0, 1))
            
        end


        function printOrderedFunctionCallList(cgt)

            A = cgt.A;

            try
                p = topological_order(A);
            catch
                fprintf('You need to install matlab BGL\n');
                return
            end

            funcCallList = {};
            for ind = 1 : numel(p)
                iprop = full(A(:, p(ind)));
                iprop = unique(iprop(iprop>0));
                if ~isempty(iprop)
                    assert(numel(iprop) == 1, 'There should be only one value for a given row');
                    propfunction = cgt.model.propertyFunctionList{iprop};
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

            % we remove duplicates (the same function can be used for several output variables)
            [~, ia, ic] = unique(funcCallList, 'first');

            ia = sort(ia); % Not sure if it is necessary, but in this way we make sure the ordering for which  funcCallList was built is respected
            funcCallList = funcCallList(ia);

            fprintf('Function call list\n');
            for ind = 1 : numel(funcCallList)
                fprintf('%s\n', funcCallList{ind});
            end
            
        end
        
            
    end
    
end
    
