classdef ImpedanceControlModel < ControlModel
%
    
    methods

        function model = ImpedanceControlModel(inputparams)

            model = model@ControlModel(inputparams);
            
        end
        
        function model = registerVarAndPropfuncNames(model)

            model = registerVarAndPropfuncNames@ControlModel(model);
            
            varnames = {};
            varnames{end + 1} = 'omega';
            
            model = model.registerVarNames(varnames);

            model = model.setAsStaticVarNames({'omega'});

            fn = @ImpedanceControlModel.setupI;
            inputvarnames = {'I'};
            outputvarname = 'controlEquation';
            model = model.registerPropFunction({outputvarname, fn, inputvarnames});
            
        end

        function state = setupI(model, state)

            state.controlEquation = state.I - 1;
            
        end

        function cleanState = addStaticVariables(model, cleanState, state)

            cleanState = model@ControlModel(model, cleanState, state);

            cleanState.omega = state.omega;

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
