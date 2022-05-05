classdef ControlModel < BaseModel

    properties
        
        Imax
        CRate
        lowerCutoffVoltage
        upperCutoffVoltage
        controlPolicy
    end
    
    
    methods

        function model = ControlModel(paramobj)
            model = model@BaseModel();
            
            fdnames = {'controlPolicy'     , ...
                       'CRate'             , ...
                       'lowerCutoffVoltage', ...
                       'upperCutoffVoltage'};
            model = dispatchParams(model, paramobj, fdnames);
            
        end
        
        function model = registerVarAndPropfuncNames(model)

            model = registerVarAndPropfuncNames@BaseModel(model);
            
            varnames = {};
            
            % Terminal voltage / [V]
            varnames{end + 1} = 'E';
            % Terminal Current / [A]
            varnames{end + 1} = 'I';
            % Control type
            varnames{end + 1} = 'ctrlType';            
            % CC_discharge
            % CC_charge
            % CV_discharge
            % CV_charge            
            
            % Terminal voltage variation/ [V/s]
            varnames{end + 1} = 'dEdt';
            % Terminal Current variation / [A/s]
            varnames{end + 1} = 'dIdt';
            
            % Equation that relates E and I (depends on geometry and discretization)
            varnames{end + 1} = 'EIequation';

            % control equation
            varnames{end + 1} = 'controlEquation';
            
        end

        function state = prepareStepControl(model, state, state0, dt, drivingForces)
        % Note : Attach to state the values necessary for the control. This is run only once at the beginning of a time step
        % default is nothing.
        end
        
        function state = updateControlEquation(model, state)
            
            Imax = model.Imax;
            Emin = model.lowerCutoffVoltage;
            Emax = model.upperCutoffVoltage;

            E = state.E;
            I = state.I;            
            ctrlVal  = state.ctrlVal;
            ctrlType = state.ctrlType;

            switch ctrlType
              case 'constantCurrent'
                ctrleq = I - ctrlVal;
              case 'constantVoltage'
                ctrleq = (E - ctrlVal)*1e5;
            end
            
            state.controlEquation = ctrleq;
                        
        end
        

        function state = updateControlAfterConvergence(model, state, state0, dt)
        % Note : This function is called in updateAfterConvergence after convergence and gives possibility to detect control switch.
        % default is nothing.
                
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
