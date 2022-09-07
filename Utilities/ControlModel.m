classdef ControlModel < BaseModel

    properties
        
        CRate
        % Control Policy (string). It can take following values
        % - 'CCCV'
        % - 'IEswitch'
        % - 'CV'
        controlPolicy
        
    end
    
    
    methods

        function model = ControlModel(paramobj)

            model = model@BaseModel();
            
            fdnames = {'controlPolicy'     , ...
                       'CRate'};
            model = dispatchParams(model, paramobj, fdnames);
            
        end
        
        function model = registerVarAndPropfuncNames(model)

            model = registerVarAndPropfuncNames@BaseModel(model);
            
            varnames = {};
            
            % Terminal voltage / [V]
            varnames{end + 1} = 'E';
            % Terminal Current / [A]
            varnames{end + 1} = 'I';
            % Equation that relates E and I (depends on geometry and discretization)
            varnames{end + 1} = 'EIequation';
            % control equation
            varnames{end + 1} = 'controlEquation';
            
            model = model.registerVarNames(varnames);

        end

        function state = prepareStepControl(model, state, state0, dt, drivingForces)
        % Note : Attach to state the values necessary for the control. This is run only once at the beginning of a time step
        % default is nothing.
        end
        
        function state = updateControlEquation(model, state)
        % Implemented by child model
        % Default do nothing
            state.controlEquation = [];
        end
        
        function state = updateControlState(model, state)
        % Implemented by child model.
        % Default do nothing (returns state unchanged).
        end
        
        function state = updateControlAfterConvergence(model, state, state0, dt)
        % Note : This function is called in updateAfterConvergence after convergence and gives possibility to detect control switch.
        % default is nothing.
                
        end

        function cleanState = addStaticVariables(model, cleanState, state)
        %  nothing is done per default
        end
        
        
    end
    
    
        
end



%{
Copyright 2021-2022 SINTEF Industry, Sustainable Energy Technology
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
