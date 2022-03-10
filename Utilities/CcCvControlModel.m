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
            varnames{end + 1} = EIequation;

            % control equation
            varnames{end + 1} = controlEquation;
            
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
