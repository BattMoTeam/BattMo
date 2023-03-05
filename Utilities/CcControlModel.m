classdef CcControlModel < ControlModel

    properties
        
        Imax
        
    end
    
    
    methods

        function model = CcControlModel(paramobj)
            
            model = model@ControlModel(paramobj);
            
        end
        
        function model = registerVarAndPropfuncNames(model)

            model = registerVarAndPropfuncNames@ControlModel(model);
            
            varnames = {};
            % control value that can be either a voltage or a current
            varnames{end + 1} = 'ctrlVal';
            % Control type (string). This value is redundant in this model with the initialControl property (see parent class ControlModel)
            % - 'charge'
            % - 'discharge'
            varnames{end + 1} = 'ctrlType';            
            model = model.registerVarNames(varnames);
            
            fn = @CcControlModel.updateControlEquation;
            model = model.registerPropFunction({'controlEquation', fn, {'ctrlVal', 'I'}});
            
        end

        
        function state = updateControlEquation(model, state)

            I        = state.I;
            ctrlVal  = state.ctrlVal;

            ctrleq = I - ctrlVal;
            
            state.controlEquation = ctrleq;
                        
        end

        function val = rampupControl(model, t, tup)

            switch model.initialControl
              case 'discharging'
                inputI = model.Imax;
              case 'charging'
                inputI = -model.Imax;
              otherwise
                error('initControl not recognized');
            end

            % hardcoded for now
            rampupcase = 'sineup';        
            switch rampupcase
                
              case 'sineup'
                val = (t <= tup) .* sineup(0, inputI, 0, tup, t) + (t > tup) .* inputI;
                
              case 'linear'
                val = (t <= tup).*t./tup .* inputI +  (t > tup) .* inputI;
                
            end
            
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
