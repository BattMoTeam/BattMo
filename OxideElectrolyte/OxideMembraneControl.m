classdef OxideMembraneControl < BaseModel
    
    properties

        controlType
        
    end
    
    
    methods
        
        function model = OxideMembraneControl(inputparams)

            model = model@BaseModel();

            fdnames = {'controlType'};
            model = dispatchParams(model, inputparams, fdnames);
            
        end

        function model = registerVarAndPropfuncNames(model)

            varnames = {};

            % Potential
            varnames{end + 1} = 'U';
            % Current
            varnames{end + 1} = 'I';            
            % Value of the variable that is controled
            varnames{end + 1} = 'ctrlVal';
            % Control Equation
            varnames{end + 1} = 'controlEquation';

            model = model.registerVarNames(varnames);
            
            fn = @OxideMembraneControl.setupControl;
            switch model.controlType
              case 'current'
                inputnames = {'I', 'ctrlVal'};
              case 'voltage'
                inputnames = {'U', 'ctrlVal'};
              otherwise
                error('controlType not recognized');
            end
            model = model.registerPropFunction({'controlEquation', fn, inputnames});

        end

        function state = setupControl(model, state)

            switch model.controlType
              case 'current'
                state.controlEquation = state.I - state.ctrlVal;
              case 'voltage'
                state.controlEquation = state.U - state.ctrlVal;
              otherwise
                error('control type not recognized');
            end
            
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
