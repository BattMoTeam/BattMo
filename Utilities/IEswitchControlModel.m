classdef IEswitchControlModel < ControlModel

    properties
        
        Imax
        lowerCutoffVoltage
        upperCutoffVoltage
        
    end
    
    
    methods

        function model = IEswitchControlModel(paramobj)
            
            model = model@ControlModel(paramobj);
            
            fdnames = {'lowerCutoffVoltage', ...
                       'upperCutoffVoltage'};
            model = dispatchParams(model, paramobj, fdnames);
            
        end
        
        function model = registerVarAndPropfuncNames(model)

            model = registerVarAndPropfuncNames@ControlModel(model);
            
            varnames = {};
            % Control type (string)
            % - 'constantCurrent'
            % - 'constantVoltage'
            varnames{end + 1} = 'ctrlType';            
            % control value that can be either a voltage or a current
            varnames{end + 1} = 'ctrlVal';            

            model = model.registerVarNames(varnames);
            
            fn = @IEswitchControlModel.updateControlEquation;
            model = model.registerPropFunction({'controlEquation', fn, {'ctrlType', 'ctrlVal', 'E', 'I'}});
            
        end

        
        function state = updateControlEquation(model, state)
            
            Imax = model.Imax;
            Emin = model.lowerCutoffVoltage;
            Emax = model.upperCutoffVoltage;

            E        = state.E;
            I        = state.I;            
            ctrlVal  = state.ctrlVal;
            ctrlType = state.ctrlType;

            switch ctrlType
              case 'constantCurrent'
                ctrleq = I - ctrlVal;
              case 'constantVoltage'
                %% TODO : fix hard-coded scaling
                ctrleq = (E - ctrlVal)*1e5;
            end
            
            state.controlEquation = ctrleq;
                        
        end
    
        function cleanState = addStaticVariables(model, cleanState, state)
            cleanState.ctrlType = state.ctrlType;
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
