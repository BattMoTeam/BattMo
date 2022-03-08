classdef ControlModel < BaseModel

    properties
        
        I0
        CRate
        lowerCutoffVoltage
        upperCutoffVoltage
        dIdtLimit
        dEdtLimit
        
    end
    
    
    methods

        function model = ControlModel(paramobj)
            model = model@BaseModel();
            
            fdnames = {'CRate'             , ...
                       'lowerCutoffVoltage', ...
                       'upperCutoffVoltage', ...
                       'dEdtLimit'         , ...
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

        function state = updateControlEquation(model, state)

            I0      = model.I0;
            Emin    = model.lowerCutoffVoltage;
            Emax    = model.upperCutoffVoltage;
            dEdtMin = model.dEdtLimit;
            dIdtMin = model.dIdtLimit;
            
            E = state.E;
            I = state.I;
            dEdt = state.dEdt;
            dIdt = state.dIdt;
            
            ctrlType = state.ctrlType;
            
            % convert to non AD variable
            Eval = value(E);
            Ival = value(I);
            
            % Note : dEdt, dIdt should be sent as non AD (so that we do not have to convert those)
            ctrlEqs.CC_discharge = I + I0;
            ctrlEqs.CC_charge    = I - I0;
            ctrlEqs.CV_discharge = I;
            ctrlEqs.CV_charge    = (E - Emax)*1e5;
               
            switch ctrlType
              case 'CC_discharge'
                if (Eval >= Emin) 
                    newCtrlType = ctrlType;
                    ctrleq = I + I0;
                elseif (Eval < Emin) 
                    newCtrlType = 'CV_discharge';
                    ctrleq = I;
                else
                    error('not detected');
                end
              case 'CC_charge'
                if (Eval <= Emax) 
                    newCtrlType = ctrlType;
                    ctrleq = I - I0;
                elseif (Eval > Emax) 
                    newCtrlType = 'CV_charge';
                    ctrleq = (E - Emax)*1e5;
                else
                    error('not detected');
                end
              case 'CV_discharge'
                if (dEdt >= dEdtMin)
                    newCtrlType = ctrlType;
                    ctrleq = I;
                elseif (dEdt < dEdtMin)
                    newCtrlType = 'CC_charge';
                    ctrleq = I - I0; 
                else
                    error('not detected');
                end    
              case 'CV_charge'
                if (dIdt <= -dIdtMin)
                    newCtrlType = ctrlType;
                    ctrleq = (E - Emax)*1e5;
                elseif (dIdt > -dIdtMin)
                    newCtrlType = 'CC_discharge';
                    ctrleq = I + I0; 
                else
                    error('not detected');
                end  
              otherwise
                error('control type ctrlType not recognized');
            end

            state.ctrlType = newCtrlType;
            state.controlEquation = ctrleq;
        end
            
        
    end
    
    
        
end
