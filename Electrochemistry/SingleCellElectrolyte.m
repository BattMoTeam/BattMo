classdef SingleCellElectrolyte < BaseModel
    
    methods
        
        function model = SingleCellElectrolyte(paramobj)
            model = model@BaseModel();
        end

        function model = registerVarAndPropfuncNames(model)

            model = registerVarAndPropfuncNames@BaseModel(model);

            varnames = {};
            % Electrolyte potential
            varnames{end + 1} = 'phi';
            % Mass conservation equation (sum of ion flux  from cathode and anode vanishes)
            varnames{end + 1} = 'massCons';

            model = model.registerVarNames(varnames);
            
        end
        
    end
    
end



%{
Copyright 2021-2023 SINTEF Industry, Sustainable Energy Technology
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
