classdef seiCcCvControlModel < ControlModel

    properties
        
        Imax
        dIdtLimit
        dEdtLimit

        % Control used initially. String that can take one of the following values
        % - 'discharging'
        % - 'charging'
        initialControl
        
    end
    
    
    methods

        function model = seiCcCvControlModel(inputparams)

            model = model@ControlModel(inputparams);
            
            fdnames = {'dEdtLimit', ...
                       'dIdtLimit', ...
                       'initialControl'};
            model = dispatchParams(model, inputparams, fdnames);
            
        end

        function model = registerVarAndPropfuncNames(model)

            model = registerVarAndPropfuncNames@ControlModel(model);
            
            varnames = {};
            % Control type : string that can take following value
            % - CC_discharge1
            % - CC_discharge2
            % - CC_charge1
            % - CV_charge2
            varnames{end + 1} = 'ctrlType';            
            varnames{end + 1} = 'time';
            model = model.registerVarNames(varnames);
            
            
            fn = @CcCvControlModel.updateControlEquation;
            model = model.registerPropFunction({'controlEquation', fn, {'E', 'I'}});
            
        end

        function state = prepareStepControl(model, state, state0, dt)
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
              case 'CC_discharge1'
                ctrleq = I - Imax;
              case 'CV_discharge2'
                ctrleq = (E - Emin);
              case 'CC_charge1'
                time = state.time;
                reftime = state.reftime;
                ctrleq = I + rampupControl(time - reftime, 1*minute, Imax);
              case 'CV_charge2'
                ctrleq = (E - Emax);
              otherwise
                error('ctrlType not recognized');
            end
            
            state.controlEquation = ctrleq;
            
        end

        function state = updateControlState(model, state, state0,dt)
            
            Emin = model.lowerCutoffVoltage;
            Imax = model.Imax;
            
            ctrlType = state.ctrlType;
            E = state.E;
            I = state.I;
            
            if strcmp(ctrlType, 'CC_discharge1')
                if E <= Emin
                    state.ctrlType = 'CV_discharge2';
                    fprintf('switch control from CC to CV\n');
                end
            end
            
            if strcmp(ctrlType, 'CV_discharge2')
                if I > Imax
                    state.ctrlType = 'CC_discharge1';
                    fprintf('switch control from CV to CC\n');
                end
            end    
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
                if (E <= Emin) 
                    nextCtrlType = 'CC_discharge2';
                end
                
              case 'CV_discharge2'
                
                nextCtrlType = 'CV_discharge2';
                if (abs(dIdt) <= dIdtMin)
                    nextCtrlType = 'CC_charge1';
                    state.reftime = state.time;
                end
                
              case 'CC_charge1'

                nextCtrlType = 'CC_charge1';
                if (E >= Emax) 
                    nextCtrlType = 'CV_charge2';
                end 
                
              case 'CV_charge2'

                nextCtrlType = 'CV_charge2';
                if (abs(dIdt) < dIdtMin)
                    nextCtrlType = 'CC_discharge1';
                end                  
                
              otherwise
                
                error('controlType not recognized');
                
            end
            
            state.nextCtrlType = nextCtrlType;
            
        end
        
        function cleanState = addStaticVariables(model, cleanState, state)
            
            cleanState = addStaticVariables@ControlModel(model, cleanState, state);
            cleanState.ctrlType = state.ctrlType;
            if isfield(state, 'reftime')
                cleanState.reftime = state.reftime;
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
