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
            varnames{end + 1} = 'prevCtrlType';
            
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
            prevCtrlType = state.prevCtrlType;
            
            % convert to non AD variable
            Eval = value(E);
            Ival = value(I);
            
            % Note : dEdt, dIdt should be sent as non AD (so that we do not have to convert those)
            ctrlEqs.CC_discharge = I + I0;
            ctrlEqs.CC_charge    = I - I0;
            ctrlEqs.CV_discharge = I;
            ctrlEqs.CV_charge    = (E - Emax)*1e5;

            switch prevCtrlType
              
              case 'CC_discharge'

                assert(ismember(ctrlType, {'CC_discharge', 'CV_discharge'}));

                if (Eval >= Emin) 
                    newCtrlType = 'CC_discharge';
                else 
                    newCtrlType = 'CV_discharge';
                end
            
              case 'CV_discharge'

                assert(ismember(ctrlType, {'CV_discharge', 'CC_charge'}));

                if (dEdt >= dEdtMin)
                    newCtrlType = 'CV_discharge';
                else
                    newCtrlType = 'CC_charge';
                end

              case 'CC_charge'
                
                assert(ismember(ctrlType, {'CC_charge', 'CV_charge'}));

                if (Eval <= Emax) 
                    newCtrlType = 'CC_charge';
                else
                    newCtrlType = 'CV_charge';
                end 
                
              case 'CV_charge'
                
                assert(ismember(ctrlType, {'CV_charge', 'CC_discharge'}));
                
                if (dIdt <= - dIdtMin)
                    newCtrlType = 'CV_charge';
                else
                    newCtrlType = 'CC_discharge';
                end                  
                
            end

            if ~strcmp(ctrlType, newCtrlType)
                fprintf('\n\n *** control switch %s -> %s\n\n', ctrlType, newCtrlType);
            end
            
            state.ctrlType = newCtrlType;
            state.controlEquation = ctrlEqs.(newCtrlType);
            
        end
            
        
    end
    
    
        
end
