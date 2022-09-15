classdef ComputationalGraphFilter

    properties
        graph
        model
        A
        nodenames
        includeNodeNames
        excludeNodeNames
    end
    
    methods
        
        function cgf = ComputationalGraphFilter(model)
            if isempty(model.propertyFunctionList)
                model = model.registerVarAndPropfuncNames();
            end
            cgf.model = model;
            g = setupGraph(model);
            cgf.graph = g;
            % In adjacency matrix A, 
            % - column index      : output variable index (as in model.varNameList)
            % - row index         : input variable (as in model.varNameList)
            % - coefficient value : property function index as in model.propertyFunctionList
            cgf.A         = adjacency(g, 'weighted');
            cgf.nodenames = g.Nodes.Variables;

            
        end

        function [propfuncs, nodenames] = findPropFunction(cgf, nodename)
            
            nodenames = cgf.nodenames;
            A         = cgf.A;
            model     = cgf.model;
            
            indSelectedNodenames = regexp(nodenames, nodename, 'once');
            indSelectedNodenames = cellfun(@(x) ~isempty(x), indSelectedNodenames);
            indPropfunctions = A(:, indSelectedNodenames);
            indPropfunctions = unique(indPropfunctions(:));

            propfuncs = model.propertyFunctionList(indPropfunctions(indPropfunctions > 0));

            if numel(propfuncs) == 1
                propfuncs = propfuncs{1};
            end
            
        end

        function printPropFunction(cgf, nodename)
            propfuncs = cgf.findPropFunction(nodename);

            if numel(propfuncs) == 1
                propfuncs = {propfuncs};
            end

            for iprop = 1 : numel(propfuncs)
                propfunc = propfuncs{iprop};
                varname = propfunc.varname;
                varname_s = varname.resolveIndex();
                for ind = 1 : numel(varname_s)
                    fullname = varname_s{ind}.getIndexedFieldname();
                    fprintf('%s : ', fullname);
                end
                fprintf('%s\n', propfuncs{iprop}.getFunctionSignature());
            end
            
            
        end
        
        function [g, edgelabels] = getFilteredGraph(cgf, varargin)
            
            opt = struct('type', 'ascendant', ...
                         'oneParentOnly', false);
            opt = merge_options(opt, varargin{:});

            propfunctionlist = cgf.model.propertyFunctionList;
            nodenames        = cgf.nodenames;
            includeNodeNames = cgf.includeNodeNames;
            excludeNodeNames = cgf.excludeNodeNames;

            A = (cgf.A);
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
        
        function [g, edgelabels] = setupAscendantGraph(cgf, varargin)
            [g, edgelabels] = cgf.setupGraph('type', 'ascendant', varargin{:});
        end
        
        function [g, edgelabels] = setupDescendantGraph(cgf, varargin)
            [g, edgelabels] = cgf.setupGraph('type', 'descendant', varargin{:});            
        end
        
    end
    
end
    
