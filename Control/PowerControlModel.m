classdef PowerControlModel < ControlModel

    properties

        dischargingPower
        chargingPower
        dischargingTime
        chargingTime
        
    end
    
    
    methods

        function model = PowerControlModel(paramobj)

            model = model@ControlModel(paramobj);
            
            fdnames = {'dischargingPower', ...
                       'chargingPower'   , ...
                       'dischargingTime' , ...
                       'chargingTime'};
            
            model = dispatchParams(model, paramobj, fdnames);
            
        end

        function model = registerVarAndPropfuncNames(model)

            model = registerVarAndPropfuncNames@ControlModel(model);
            
            varnames = {};
            % Control type : string that can take following value
            % - discharge
            % - charge
            varnames{end + 1} = 'ctrlType';
            varnames{end + 1} = 'time';
            model = model.registerVarNames(varnames);
            
            fn = @PowerControlModel.updateControlEquation;
            model = model.registerPropFunction({'controlEquation', fn, {'ctrlType', 'E', 'I'}});
            
        end


        function state = updateControlEquation(model, state)
            
            Pcharge    = model.chargingPower;
            Pdischarge = model.dischargingPower;
            
            E = state.E;
            I = state.I;            
            ctrlType = state.ctrlType;
            
            switch ctrlType
              case 'discharge'
                ctrleq = E*I - Pdischarge;
              case 'charge'
                ctrleq = E*I + Pcharge;
              otherwise
                error('ctrlType not recognized');
            end
            state.controlEquation = ctrleq;
            
        end

        function state = prepareStepControl(model, state, state0, dt, drivingForces)

            cT = model.chargingTime;
            dT = model.dischargingTime;
            
            time = state.time;

            r = time - floor(time/(dT + cT))*(dT + cT);
            
            if r < cT
                state.ctrlType = 'charge';
            else
                state.ctrlType = 'discharge';
            end
            
        end

        function cleanState = addStaticVariables(model, cleanState, state)
            cleanState.ctrlType = state.ctrlType;
        end
        
    end
    
    
        
end



%{
Copyright 2021-2023 SINTEF Industry, Sustainable Energy Technology
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
