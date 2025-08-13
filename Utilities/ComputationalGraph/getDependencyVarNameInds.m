function [depvarnameinds, propfuncinds, propvarnameinds, propdeplevels, rootinds, rootdeplevels] = getDependencyVarNameInds(varnameinds, A)
% The inputs are
% - varnameinds : an array of indices (indices must  A)
% - A           : a square matrix corresponding to the adjency matrix.
%                 The values of of the coefficient in the matrix corresponds to the index of the PropFunction, see ComputationalGraph for details. 
%
% The outputs are
% - depvarnameinds  : Indices of the variables that are required for the evaluation of the variables given by varnameinds
%                     for the functional dependencies given by the adjency matrix A.
%                     It is obtain by concatenating rootinds and propvarnameinds
% - propfuncinds    : Indices of PropFunctions that will be evaluated to obtain the depvarnameinds
% - propvarnameinds : Indices of the variables that corresponds to the PropFunction indices given in propfuncinds,
%                     Same size as propfuncinds
% - propdeplevels   : Level of dependency (increases by one for each evaluation of a PropFunction that is needed,
%                     before requiring the evaluation of the current PropFunction),
%                     Same size as propfuncinds
% - rootinds        : Indices that depends on the indices given by varnameinds throught the adjency matrix A, but which corresponds to root variable
% - rootdeplevels   : Level of dependency (increases by one for each evaluation of a PropFunction that is needed to obtain the PropFunction given by propfuncinds)
%                     Same size as rootinds
% 
    [propfuncinds   , ...
     propvarnameinds, ...
     propdeplevels  , ...
     rootinds     , ...
     rootdeplevels] =  getDependencyVarNameIndsRecursive(varnameinds, ...
                                                           []       , ... 
                                                           []       , ... 
                                                           []       , ... 
                                                           []       , ... 
                                                           []       , ...
                                                           0        , ... 
                                                           A);

    depvarnameinds = [rootinds; propvarnameinds];
    
end

function [propfuncinds   , ...
          propvarnameinds, ...
          propdeplevels  , ...
          rootinds       , ...
          rootdeplevels] = getDependencyVarNameIndsRecursive(varnameinds      , ...
                                                               propfuncinds   , ...
                                                               propvarnameinds, ...
                                                               propdeplevels  , ...
                                                               rootinds       , ...
                                                               rootdeplevels  , ...
                                                               level          , ...
                                                               A)
% Helper function for recursive computation

    depvarnameinds = [];
    
    for ivar = 1 : numel(varnameinds)

        varnameind = varnameinds(ivar);
        c = A(:, varnameind);

        if all(c == 0)
            %  no dependency
            rootinds      = vertcat(rootinds, varnameind);
            rootdeplevels = vertcat(rootdeplevels, level);
            
        else
            % the non-zero entry in the column is is unique and give the index of the property function
            
            propfuncind = max(c);

            depvarnameinds  = vertcat(depvarnameinds, find(c));
            propfuncinds    = vertcat(propfuncind, propfuncinds);
            propvarnameinds = vertcat(varnameind, propvarnameinds);
            propdeplevels   = vertcat(level, propdeplevels);
            
        end
        
    end        
    
    if ~isempty(depvarnameinds)

        varnameinds = depvarnameinds;
        
        [propfuncinds   , ...
         propvarnameinds, ...
         propdeplevels  , ...
         rootinds     , ...
         rootdeplevels] = getDependencyVarNameIndsRecursive(varnameinds    , ...
                                                            propfuncinds   , ...
                                                            propvarnameinds, ...
                                                            propdeplevels  , ...
                                                            rootinds       , ...
                                                            rootdeplevels  , ...
                                                            level + 1      , ...
                                                            A); 
        
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
