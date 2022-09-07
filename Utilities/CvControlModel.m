classdef CvControlModel < ControlModel
    
    properties
        inputVoltage
    end
    
    methods


        function model = CvControlModel(paramobj)
            
            model = model@ControlModel(paramobj);
            
            fdnames = {'inputVoltage'};
            model = dispatchParams(model, paramobj, fdnames);
            
        end
        
        function model = registerVarAndPropfuncNames(model)

            fn = @CvControlModel.updateControlEquation;
            model = model.registerPropFunction({'controlEquation', fn, {'E'}});
            
        end
        
        function state = updateControlEquation(model, state)
            
            inputE  = model.inputVoltage;

            E = state.E;

            %% TODO : fix scaling
            ctrleq = (E - inputE)*1e5;
            
            state.controlEquation = ctrleq;
                        
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
