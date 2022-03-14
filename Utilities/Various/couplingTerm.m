classdef couplingTerm
% Structure that describes the topology of a coupling: Two components (each with their own grid structure) are coupled
% through cells or faces
%
% The nature of the coupling (which physical process it corresponds to) is not addressed here.
%
    properties
        name           % The coupling is given an name (string)
        
        componentnames % List of the component names that enter in this coupling (cell array of string)
        
        couplingcells  % Array of cells : For everey i, we have a coupling from the cell indexed by couplingcells(i, 1) in
                       % component 1 with cell indexed by couplingcells(i, 2) in the component 2. The ordering of the
                       % component follows the ordering in componentnames
        
        couplingfaces  % Same as couplingcells but for faces (need not be assigned - it depends on how the coupling is used in
                       % practice)
    end
    
    methods
        
        function obj = couplingTerm(name, compnames)
            obj.name = name;
            obj.componentnames = compnames;
        end
        
    end
    
end



%{
Copyright 2021-2022 SINTEF Industry, Sustainable Energy Technology
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
