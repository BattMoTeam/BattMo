classdef ComputationalGraphFilter

    properties
        graph
        A
        nodenames
        includeNodeNames
        excludeNodeNames
    end
    
    methods
        
        function cgf = ComputationalGraphFilter(graph)
            cgf.graph     = graph;
            cgf.A         = adjacency(graph);
            cgf.nodenames = graph.Nodes.Variables;
        end

        
        function g = setupGraph(cgf, varargin)
            
            opt = struct('type', 'ascendant');
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
            
            nodes = getNodeDependencyListByName(includeNodeNames, nodenames, A);

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
        
        function g = setupAscendantGraph(cgf)
            g = cgf.setupGraph('type', 'ascendant');
        end
        
        function g = setupDescendantGraph(cgf)
            g = cgf.setupGraph('type', 'descendant');            
        end
        
    end
    
end
    
