classdef ConcentrationReactionThermalModel < BaseModel

    properties

        Masses
        Thermal
        
    end
    
    methods

        function model = ConcentrationReactionThermalModel()

            model.Masses  = ConcentrationReactionModel();
            model.Thermal = ThermalModel();
            
        end

        function model = registerVarAndPropfuncNames(model)
            
            %% Declaration of the Dynamical Variables and Function of the model

            model = registerVarAndPropfuncNames@BaseModel(model);

            fn = @ConcentrationReactionThermalModel.updateOCP;
            inputnames = {{'Masses', 'Reaction', 'c_s'}, {'Thermal', 'T'}};
            model = model.registerPropFunction({{'Masses', 'Reaction', 'OCP'}, fn, inputnames});

            fn = @ConcentrationReactionThermalModel.updateThermalSource;
            inputnames = {{'Masses', 'Reaction', 'R'}};
            model = model.registerPropFunction({{'Thermal', 'source'}, fn, inputnames});
            
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
