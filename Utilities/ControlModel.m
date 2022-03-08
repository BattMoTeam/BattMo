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
            % CC-discharge
            % CC-charge
            % CV-discharge
            % CV-charge            
            
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
            dVdtMin = model.dEdtLimit;
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
            
            switch ctrlType
              case 'CC-discharge'
                if (Ival == -I0) & (E >= Emin) 
                    newCtrlType = ctrlType;
                    ctrleq = I + I0;
                elseif (Ival == -I0) & (E < Emin) 
                    newCtrlType = 'CV-discharge';
                    ctrleq = I;
                end
              case 'CC-charge'
                if (Ival == I0) & (E <= Emax) 
                    newCtrlType = ctrlType;
                    ctrleq = I - I0;
                elseif (Ival == I0) & (E > Emax) 
                    newCtrlType = 'CV-charge';
                    ctrleq = (E - Emax)*1e5;
                end
              case 'CV-discharge'
                if (I == 0) & (dVdt >= dVdtMin)
                    newCtrlType = crtlType;
                    ctrleq = I;
                else (I == 0) & (dIdt < dVdtMin)
                    newCtrlType = 'CC-charge';
                    ctrleq = I - I0; 
                end    
              case 'CV-charge'
                if (Eval == Emax) & (dIdt >= dIdtMin)
                    newCtrlType = crtlType;
                    ctrleq = (E - Emax)*1e5;
                else (Eval == Emax) & (dIdt < dIdtMin)
                    newCtrlType = 'CC-charge';
                    ctrleq = I - I0; 
                end                    
            end
            
            state.ctrlType = newCtrlType;
            state.controlEquation = ctrleq;
        end
            
        
    end
    
    
        
end
