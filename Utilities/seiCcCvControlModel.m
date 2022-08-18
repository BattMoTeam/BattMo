classdef seiCcCvControlModel < CcCvControlModel


    methods
        
        function state = updateControlEquation(model, state)
            
            Imax = model.Imax;
            Emin = model.lowerCutoffVoltage;
            Emax = model.upperCutoffVoltage;
            
            E = state.E;
            I = state.I;            
            ctrlType = state.ctrlType;
            
            switch ctrlType
              case 'CC_discharge1'
                ctrleq = I - Imax;
              case 'CV_discharge2'
                ctrleq = (E - Emin);
              case 'CC_charge1'
                    ctrleq = I + 1e-2*Imax;
              case 'CC_charge2'
                    ctrleq = I + Imax;
              case 'CV_charge3'
                ctrleq = (E - Emax);
              otherwise
                error('ctrlType not recognized');
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
              
              case 'CC_discharge1'
                
                nextCtrlType = 'CC_discharge1';
                if (E <= 1.3*Emin) 
                    nextCtrlType = 'CV_discharge2';
                end
            
              case 'CV_discharge2'
                
                nextCtrlType = 'CV_discharge2';
                if (abs(dIdt) <= dIdtMin)
                    nextCtrlType = 'CC_charge1';
                end
              
              case 'CC_charge1'

                nextCtrlType = 'CC_charge1';
                if (E >= 1.3*Emin) 
                    nextCtrlType = 'CC_charge2';
                end 
            
              case 'CC_charge2'

                nextCtrlType = 'CC_charge2';
                if (E >= Emax) 
                    nextCtrlType = 'CV_charge3';
                end 
                
              case 'CV_charge3'

                nextCtrlType = 'CV_charge3';
                if (abs(dIdt) < dIdtMin)
                    nextCtrlType = 'CC_discharge1';
                end                  
                
              otherwise
                
                error('controlType not recognized');
                
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
