function p = compatible_topological_order(A)
%
% We implement here topological sorting for a graph using the Kahn's algorithm (https://en.wikipedia.org/wiki/Topological_sorting#Kahn's_algorithm)
%
% An alternative is to use topological_order from matlab_blg package but it is not compatible with Octave.
%     
% The input A denotes the adjacency matrix of the graph
%   
% In adjacency matrix A, 
% - column index      : target variable
% - row index         : source variable
% - coefficient value : weight (not used here)
    
    nnodes = size(A, 2);
    indegree = zeros(nnodes, 1);

    for inode = 1 : nnodes

        neighbors = find(A(inode, :) > 0);
        indegree(neighbors) = indegree(neighbors) + 1;

    end

    orderedNodes = [];

    done     = false;
    iscyclic = false;

    while ~done

        rootnodes = find(indegree == 0);

        if isempty(rootnodes)
            
            iscyclic = true;
            done = true;

        else


            indegree(rootnodes) = NaN;

            
            for iroot = 1 : numel(rootnodes)

                rootnode = rootnodes(iroot);
                neighbors = find(A(rootnode, :) > 0);
                indegree(neighbors) = indegree(neighbors) - 1;
                
            end
            
            orderedNodes = vertcat(orderedNodes, rootnodes);

            if all(isnan(indegree))
                done = true;
            end

        end
        
    end

    if iscyclic
        p = [];
    else
        p = orderedNodes;
    end

end



%{
Copyright 2021-2024 SINTEF Industry, Sustainable Energy Technology
and SINTEF Digital, Mathematics & Cybernetics.

This file is part of The Battery Modeling Toolbox BattMo

BattMo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BattMo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BattMo.  If not, see <http://www.gnu.org/licenses/>.
%}
