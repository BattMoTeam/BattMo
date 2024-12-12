classdef GenericControlModel < ControlModel
  %
  properties
    
    controlsteps % cell array, each cell describing a control (see GenericControlModel.schema.json)
    
  end
  
  methods
    
    function model = GenericControlModel(inputparams)
      
      model = model@ControlModel(inputparams);
      
      fdnames = {'controlsteps'};
      model = dispatchParams(model, inputparams, fdnames);
      
    end
    
    function model = registerVarAndPropfuncNames(model)
      
      model = registerVarAndPropfuncNames@ControlModel(model);
      
      varnames = {};
      % Control type : string that can take following value
      % - 'rest'
      % - 'current'
      % - 'voltage'
      varnames{end + 1} = 'ctrlType';
      varnames{end + 1} = 'ctrlStepNumber';
      
      model = model.registerVarNames(varnames);
      
      model = model.setAsStaticVarName('ctrlStepNumber');
      
      % Register the functions
      fn = @CcCvControlModel.updateControlType;
      model = model.registerPropFunction({'ctrlType', fn, {'ctrlStepNumber'}});
      
      fn = @CcCvControlModel.updateControlEquation;
      model = model.registerPropFunction({'controlEquation', fn, {'ctrlType', 'E', 'I'}});
      
    end
    
    function cleanState = addStaticVariables(model, cleanState, state)
      
      cleanState.ctrlStepNumber = state.ctrlStepNumber;
      
    end
    
    function func = setupControlFunction(model)
      
      func = [];
      
    end
    
    function control = setupScheduleControl(model)
      
      control = setupScheduleControl@ControlModel(model);
      control.Generic = true;
      
    end
    
    function step = setupScheduleStep(model, timeSteppingParams)
      
      % Setup and a return the step structure that is part of the schedule which is used as input for
      % :mrst:`simulateScheduleAD`. For some control type, there is a natural construction for this structure. This is
      % why we include this method here, for convenience. It can be overloaded by derived classes. The
      % timeSteppingParams structure by default is given by the data described in :battmofile:`Utilities/JsonSchemas/TimeStepping.schema.json`
      
      step.val     = 10*minute*ones(600, 1);
      step.control = ones(600, 1);
      
    end
    
    function state = updateControlType(model, state)
      
      istep = state.ctrlStepNumber;
      
      ctrlstep = model.controlsteps{istep};
      
      state.ctrlType = ctrlstep.controltype;
      
    end
    
    function state = updateControlState(model, state, state0, dt)
      
      % This function is called by the solver at the end of a Newton step.  We check if we have triggered the
      % "termination" criteria from the control If the termination criteria is met, we increment the control step by
      % one, and proceed with the Newton algorithm
      
      ctrlType  = state.ctrlType;
      ictrlstep = state.ctrlStepNumber;
      
      ctrlstep = model.controlsteps{ictrlstep};
      termination = ctrlstep.termination;
      
      doswitch = false;
      
      switch termination.terminationtype
        case 'current'
          I = state.I;
          tI = termination.value;
          if abs(I) < tI
            doswitch = true;
          end
        case 'voltage'
          E = state.E;
          switch ctrlstep.direction
            case 'discharge'
              Emin = termination.value;
              if E < Emin
                doswitch = true;
              end
            case 'charging'
              Emax = termination.value
              if E > Emax
                doswitch = true;
              end
            otherwise
              error('direction not recognized');
          end
        case 'time'
        otherwise
          error('termination type not recognized')
      end
      
      if doswitch
        state.ctrlStepNumber = state.ctrlStepNumber + 1;
        state = model.updateValueFromControl(state);
      end
      
    end
    
    function  [arefulfilled, state] = checkConstraints(model, state, state0, dt)
      
      % This function is called when the Newton method has converged, but we want to check if the solution we have
      % obtained does not break the termination criteria.
      
    end
    
    
    function state = updateValueFromControl(model, state)
      % From the given control type, set the corresponding control variable to the expected value
      
      ictrlstep = state.ctrlStepNumber;
      controlstep = model.controlsteps{ictrlstep};
      ctrlType = controlstep.controltype;
      
      state.ctrlType = ctrlType;
      
      switch ctrlType
        
        case 'rest'
          
          state.I = 0
          
        case 'current'
          
          givenI = controlstep.value;
          state.I = givenI;
          
        case 'voltage'
          
          givenE = controlstep.value;
          state.E = givenE;
          
        otherwise
          
          error('control type not recognized');
          
      end
      
    end
    
    function state = updateControlEquation(model, state)
        
      ctrlType  = state.ctrlType;
      ictrlstep = state.ctrlStepNumber;
      
      controlstep = model.controlsteps{ictrlstep};
      
      switch ctrlType
        
        case 'rest'
          
          ctrleq = state.I
          
        case 'current'
          
          givenI = controlstep.value;
          
          ctrleq = state.I - givenI;
          
        case 'voltage'
          
          givenE = controlstep.value;
          
          ctrleq = (state.E - givenE)*1e5;
          
        otherwise
          
          error('control type not recognized');
          
      end
      
      state.controlEquation = ctrleq;
      
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
