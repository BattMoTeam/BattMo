function model = convertModelGrids(model)
% Convert the grids from sub-grid representation to a concrete grid with the data needed for a julia simulation.
    
    if isprop(model, 'G') && ~isempty(model.G)
        G = model.G.mrstFormat;
        op.haltTrans = 
        G.operators = op;
        model.G = G;
    end

    submodelnames = model.getSubModelNames();

    if ~isempty(submodelnames)

        for isub = 1 : numel(submodelnames)

            submodelname = submodelnames{isub};

            model.(submodelname) = convertModelGrids(model.(submodelname));

        end

    end
    
end

%{
Copyright 2021-2026 SINTEF Industry, Sustainable Energy Technology
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
