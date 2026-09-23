classdef ReactionModel < BaseModel

    methods

        function model = registerVarAndPropfuncNames(model)
            
            model = registerVarAndPropfuncNames@BaseModel(model);
            
            varnames = {};
            % potential in electrode
            varnames{end + 1} = 'phi_s';
            % charge carrier concentration in electrode - value at surface
            varnames{end + 1} = 'c_s';
            % potential in electrolyte
            varnames{end + 1} = 'phi_e';
            % charge carrier concentration in electrolyte
            varnames{end + 1} = 'c_e';
            % Electrode over potential
            varnames{end + 1} = 'eta';
            % Intercalation flux [mol s^-1 m^-2]
            varnames{end + 1} = 'R';
            % OCP [V]
            varnames{end + 1} = 'OCP';
            % Reaction rate coefficient [A m^-2]
            varnames{end + 1} = 'j';
            
            model = model.registerVarNames(varnames);
            
            fn = @ReactionModel.updateReactionRateCoefficient;
            inputnames = {'c_e', 'c_s'};
            model = model.registerPropFunction({'j', fn, inputnames});

            fn = @ReactionModel.updateOCP;
            inputnames = {'c_s'};
            model = model.registerPropFunction({'OCP', fn, inputnames});
            
            fn = @ReactionModel.updateEta;
            inputnames = {'phi_e', 'phi_s', 'OCP'};            
            model = model.registerPropFunction({'eta', fn, inputnames});

            fn = @ReactionModel.updateReactionRate;
            inputnames = {'eta', 'j'};
            model = model.registerPropFunction({'R', fn, inputnames});
            
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
