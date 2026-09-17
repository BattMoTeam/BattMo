classdef ConcentrationModel < BaseModel

    methods

        function model = registerVarAndPropfuncNames(model)
            
            %% Declaration of the Dynamical Variables and Function of the model
            % (setup of varnameList and propertyFunctionList)

            model = registerVarAndPropfuncNames@BaseModel(model);
            
            varnames = {};
            % Concentration
            varnames{end + 1} = 'c';
            % Mass accumulation term
            varnames{end + 1} = 'massAccum';
            % Source term
            varnames{end + 1} = 'source';
            % Mass conservation equation
            varnames{end + 1} = 'massCons';
            
            model = model.registerVarNames(varnames);
            
            fn = @ConcentrationModel.updateMassAccum;
            inputnames = {'c'};
            model = model.registerPropFunction({'massAccum', fn, inputnames});

            fn = @ConcentrationModel.updateMassCons;
            inputnames = {'massAccum', 'source'};
            model = model.registerPropFunction({'massCons', fn, inputnames});
            
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
