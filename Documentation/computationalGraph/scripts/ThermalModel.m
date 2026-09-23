classdef ThermalModel < BaseModel

    methods

        function model = registerVarAndPropfuncNames(model)
            
            model = registerVarAndPropfuncNames@BaseModel(model);
            
            varnames = {};
            % Temperature
            varnames{end + 1} = 'T';
            % Accumulation term
            varnames{end + 1} = 'accumTerm';
            % Flux flux
            varnames{end + 1} = 'flux';
            % Heat source
            varnames{end + 1} = 'source';
            % Energy conservation
            varnames{end + 1} = 'energyCons';
            
            model = model.registerVarNames(varnames);
            
            fn = @ThermalModel.updateFlux;
            model = model.registerPropFunction({'flux', fn, {'T'}});

            fn = @ThermalModel.updateAccumTerm;
            model = model.registerPropFunction({'accumTerm', fn, {'T'}});
            
            fn = @ThermalModel.updateEnergyCons;
            inputnames = {'accumTerm', 'flux', 'source'};            
            model = model.registerPropFunction({'energyCons', fn, inputnames});

            
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
