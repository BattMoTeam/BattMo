classdef CcCvControlModel < ControlModel

    properties
        
        dIdtLimit
        dEdtLimit
        
    end
    
    
    methods

        function model = CcCvControlModel(paramobj)

            model = model@ControlModel(paramobj);
            
            fdnames = {'dEdtLimit'         , ...
                       'dIdtLimit'};
            model = dispatchParams(model, paramobj, fdnames);
            
        end

        function model = registerVarAndPropfuncNames(model)

            model = registerVarAndPropfuncNames@BaseModel(model);
            
            fn = @CcCvControlModel.updateControlEquation;
            model = model.registerPropFunction({'controlEquation', fn, {'E', 'I'}});
            
        end

        function state = prepareStepControl(model, state, state0, dt, drivingForces)
            state.ctrlType = state0.nextCtrlType;
        end
        
        function state = updateControlEquation(model, state)
            
            Imax = model.Imax;
            Emin = model.lowerCutoffVoltage;
            Emax = model.upperCutoffVoltage;
            
            E = state.E;
            I = state.I;            
            ctrlType = state.ctrlType;
            
            switch ctrlType
              case 'CC_discharge'
                ctrleq = I - Imax;
              case 'CV_discharge'
                ctrleq = I;
              case 'CC_charge'
                if (value(E) <= Emax)
                    ctrleq = I + Imax;
                else
                    ctrleq = (E - Emax)*1e5;
                    state.ctrlType = 'CV_charge';
                end
              case 'CV_charge'
                ctrleq = (E - Emax)*1e5;
            end
            
            state.controlEquation = ctrleq;
            
        end
        

        function state = updateControlAfterConvergence(model, state, state0, dt)

            Imax    = model.Imax;
            Emin    = model.lowerCutoffVoltage;
            Emax    = model.upperCutoffVoltage;
            dEdtMin = model.dEdtLimit;
            dIdtMin = model.dIdtLimit;
            
            E = state.E;
            
            dEdt = (state.E - state0.E)/dt;
            dIdt = (state.I - state0.I)/dt;

            ctrlType = state.ctrlType;
            
            switch ctrlType
              
              case 'CC_discharge'
                
                if (E >= Emin) 
                    nextCtrlType = 'CC_discharge';
                else 
                    nextCtrlType = 'CV_discharge';
                end
            
              case 'CV_discharge'

                if (dEdt >= dEdtMin)
                    nextCtrlType = 'CV_discharge';
                else
                    nextCtrlType = 'CC_charge';
                end

              case 'CC_charge'

                if (E <= Emax) 
                    nextCtrlType = 'CC_charge';
                else
                    nextCtrlType = 'CV_charge';
                end 
                
              case 'CV_charge'
                
                if (dIdt >= dIdtMin)
                    nextCtrlType = 'CV_charge';
                else
                    nextCtrlType = 'CC_discharge';
                end                  
                
            end
            
            state.nextCtrlType = nextCtrlType;
            
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
