classdef GenericControlModel < ControlModel
%
    properties

        controlsteps
        
    end
    
    methods

        function model = GenericControlModel(inputparams)

            model = model@ControlModel();
            
        end

        function model = registerVarAndPropfuncNames(model)

            model = registerVarAndPropfuncNames@ControlModel(model);
            
            varnames = {};
            % Control type : string that can take following value
            % - 'rest'
            % - 'current'
            % - 'voltage'
            varnames{end + 1} = 'ctrlType';
            varnames{end + 1} = 'ctrlStepNumber';
            
            model = model.registerVarNames(varnames);

            % Register the functions
            fn = @CcCvControlModel.updateControlType;
            model = model.registerPropFunction({'ctrlType', fn, {'ctrlStepNumber'}});
            
            fn = @CcCvControlModel.updateControlEquation;
            model = model.registerPropFunction({'controlEquation', fn, {'ctrlType', 'E', 'I'}});
            
        end

        function state = updateControlType(model, state)

            ctrlstep = model.controlsteps{istep}

            istep = state.ctrlStepNumber;
            
            state.ctrlType = ctrlstep.controlstep;
            
        end

        
        function state = updateControlState(model, state, state0, dt)

            ctrlType = state.ctrlType;

            switch ctrlType

              case 'rest'

              case 'current'

              case 'voltage'

              otherwise
                
                error('control type not recognized');
                
            end
            
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
