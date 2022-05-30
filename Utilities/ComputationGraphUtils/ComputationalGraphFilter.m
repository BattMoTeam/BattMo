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

        
        function g = setupGraph(cgf)
            
            nodenames        = cgf.nodenames;
            A                = cgf.A;
            includeNodeNames = cgf.includeNodeNames;
            excludeNodeNames = cgf.excludeNodeNames;
            
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
            
            g = digraph(A, nodenames);
            
        end
        
    end
    
end
    
