classdef ConcentrationReactionModel < BaseModel

    properties

        Reaction
        Solid
        Elyte
        
    end
    
    methods

        function model = ConcentrationReactionModel()

            model.Reaction = ReactionModel();
            model.Solid    = ConcentrationModel();
            model.Elyte    = ConcentrationModel();
            
        end
        
        function model = registerVarAndPropfuncNames(model)
            
            %% Declaration of the Dynamical Variables and Function of the model

            model = registerVarAndPropfuncNames@BaseModel(model);

            fn = @ConcentrationReactionModel.updateConcentrationSource;
            inputnames = {{'Reaction', 'R'}};
            model = model.registerPropFunction({{'Solid', 'source'}, fn, inputnames});
            model = model.registerPropFunction({{'Elyte', 'source'}, fn, inputnames});
            
            fn = @ConcentrationReactionModel.updateReactionConcentrationE;
            inputnames = {{'Elyte', 'c'}};
            model = model.registerPropFunction({{'Reaction', 'c_e'}, fn, inputnames});

            fn = @ConcentrationReactionModel.updateReactionConcentrationS;
            inputnames = {{'Solid', 'c'}};
            model = model.registerPropFunction({{'Reaction', 'c_s'}, fn, inputnames});

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
