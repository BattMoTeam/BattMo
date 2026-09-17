function model = convertModelForJuliaBridge(model)

    % convert the grids from sub-grids, adding the transmissibilities and other quantities needed in the julia
    % simulation.
    model = convertModelGrids(model);

    model = convertModelForJuliaBridge_(model);
    
end


function model = convertModelForJuliaBridge_(model)

    props = properties(model);

    for iprop = 1 : numel(props)
        prop = props{iprop};
        if isa(model.(prop), 'function_handle')
            model.(prop) = func2str(model.(prop));
        end
    end

    submodelnames = model.getSubModelNames();

    if ~isempty(submodelnames)

        for isub = 1 : numel(submodelnames)

            submodelname = submodelnames{isub};
            model.(submodelname) = convertModelForJuliaBridge_(model.(submodelname));

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
