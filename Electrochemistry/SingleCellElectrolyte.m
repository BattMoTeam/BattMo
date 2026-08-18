classdef SingleCellElectrolyte < BaseModel


    properties

        cInit % concentration in electrolyte
        phiInit % concentration in electrolyte

        volumeFraction
        
    end
    
    methods
        
        function model = SingleCellElectrolyte(inputparams)
            model = model@BaseModel();

            model.G     = inputparams.G;
            model.cInit = inputparams.cInit;
            model.phiInit = inputparams.phiInit;
        end

        function model = registerVarAndPropfuncNames(model)

            model = registerVarAndPropfuncNames@BaseModel(model);

            varnames = {};
            % Electrolyte potential
            varnames{end + 1} = 'phi';
            % Electrolyte potential
            varnames{end + 1} = 'c';
            % Mass conservation equation (sum of ion flux  from cathode and anode vanishes)
            varnames{end + 1} = 'massCons';

            model = model.registerVarNames(varnames);

            fn = @SingleCellElectrolyte.updateConcentration;
            model = model.registerPropFunction({'c', fn, {}});
            
            fn = @SingleCellElectrolyte.updatePotential;
            model = model.registerPropFunction({'phi', fn, {}});            
        end

        function state = updateConcentration(model, state)

            state.c = model.cInit;
        end
        
        function state = updatePotential(model, state)

            state.phi = model.phiInit;
        end
        
        
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
