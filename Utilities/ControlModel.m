classdef ControlModel < BaseModel

    properties
        
        Imax
        CRate
        lowerCutoffVoltage
        upperCutoffVoltage
        controlPolicy
    end
    
    
    methods

        function model = ControlModel(paramobj)
            model = model@BaseModel();
            
            fdnames = {'controlPolicy'     , ...
                       'CRate'             , ...
                       'lowerCutoffVoltage', ...
                       'upperCutoffVoltage'};
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
        % Note : Attach to state the values necessary for the control. This is run only once at the beginning of a time step
        % default is nothing.
        end
        
        function state = updateControlEquation(model, state)
            
            Imax = model.Imax;
            Emin = model.lowerCutoffVoltage;
            Emax = model.upperCutoffVoltage;

            E = state.E;
            I = state.I;            
            ctrlval  = state.ctrlval;
            ctrltype = state.ctrltype;

            switch ctrltype
              case 'I'
                eqs{end + 1} = I - ctrlval;
              case 'E'
                eqs{end + 1} = (E - ctrlval)*1e5;
            end
            
            state.controlEquation = ctrleq;
                        
        end
        

        function state = updateControlAfterConvergence(model, state, state0, dt)
        % Note : This function is called in updateAfterConvergence after convergence and gives possibility to detect control switch.
        % default is nothing.
                
        end
            
        
    end
    
    
        
end
