classdef GenericExampleModel < BaseModel

    methods
        
        function model = registerVarAndPropfuncNames(model)
            
            model = registerVarAndPropfuncNames@BaseModel(model);

            n = 3;
            p = 2;
            q = 3;
            
            varnames = {};
            varnames{end + 1} = VarName({}, 'x', n);
            varnames{end + 1} = VarName({}, 'y', p);
            varnames{end + 1} = VarName({}, 'f', q);
            
            model = model.registerVarNames(varnames);
            
            fn = @ThermalModel.updateY1;
            outputvarname = VarName({}, 'y', p, 1);
            inputvarnames = {VarName({}, 'x', n, 1 : 2), VarName({}, 'y', n, 2)};

            model = model.registerPropFunction({outputvarname, fn, inputvarnames});

            fn = @ThermalModel.updateY2;
            outputvarname = VarName({}, 'y', p, 2);
            inputvarnames = {VarName({}, 'x', n, 1 : 3)};
            model = model.registerPropFunction({outputvarname, fn, inputvarnames});

            fn = @ThermalModel.updateF1;
            outputvarname = VarName({}, 'f', q, 1);
            inputvarnames = {VarName({}, 'y', p, 1 : 2)};
            model = model.registerPropFunction({outputvarname, fn, inputvarnames});

            fn = @ThermalModel.updateF2;
            outputvarname = VarName({}, 'f', q, 2);
            inputvarnames = {VarName({}, 'y', p, 1 : 2), VarName({}, 'x', n, 2)};
            model = model.registerPropFunction({outputvarname, fn, inputvarnames});

            fn = @ThermalModel.updateF3;
            outputvarname = VarName({}, 'f', q, 3);
            inputvarnames = {VarName({}, 'x', p, [1, 3])};
            model = model.registerPropFunction({outputvarname, fn, inputvarnames});

            
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
