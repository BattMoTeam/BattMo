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
            cgf.model       = model;
            [g, edgelabels] = setupGraph(model);
            cgf.graph       = g;
            cgf.A           = adjacency(g, 'weighted');
            cgf.nodenames   = g.Nodes.Variables;
        end


        function propfuncs = findPropFunction(cgf, nodename)
            
            nodenames = cgf.nodenames;
            A         = cgf.A;
            model     = cgf.model;
            
            indSelectedNodenames = regexp(nodenames, nodename, 'once');
            indSelectedNodenames = cellfun(@(x) ~isempty(x), indSelectedNodenames);
            indPropfunctions = A(:, indSelectedNodenames);
            indPropfunctions = unique(indPropfunctions(:));

            propfuncs = model.propertyFunctionList(indPropfunctions(indPropfunctions > 0));
            
        end
        
        function g = setupGraph(cgf, varargin)
            
            opt = struct('type', 'ascendant', ...
                         'oneParentOnly', false);
            opt = merge_options(opt, varargin{:});
            
            nodenames        = cgf.nodenames;
            includeNodeNames = cgf.includeNodeNames;
            excludeNodeNames = cgf.excludeNodeNames;

            A = (cgf.A);
            if strcmp(opt.type, 'descendant')
                A = A';
            end
                        
            
            if isempty(includeNodeNames)
                g = cgf.graph;
                return
            end
            
            nodes = getNodeDependencyListByName(includeNodeNames, nodenames, A, 'oneParentOnly', opt.oneParentOnly);

            nodenames = nodenames(nodes);
            A = A(nodes, nodes);
            
            if ~isempty(excludeNodeNames)
                removenodes = regexp(nodenames, excludeNodeNames, 'once');
                removenodes = cellfun(@(x) ~isempty(x), removenodes);
                nodenames = nodenames(~removenodes);
                A = A(~removenodes, ~removenodes);
            end
            
            if strcmp(opt.type, 'descendant')
                A = A';
            end
            
            g = digraph(A, nodenames);
            
        end
        
        function g = setupAscendantGraph(cgf, varargin)
            g = cgf.setupGraph('type', 'ascendant', varargin{:});
        end
        
        function g = setupDescendantGraph(cgf, varargin)
            g = cgf.setupGraph('type', 'descendant', varargin{:});            
        end
        
    end
    
end
    
