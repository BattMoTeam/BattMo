classdef ControlModel < BaseModel
%
% Base class for a control model
%
%
    properties

        %
        % C Rate
        %
        CRate

        %
        % Control policy (string). It can take following values
        %
        % - 'CCDischarge'
        % - 'CCCharge'
        % - 'CC'
        % - 'CCCV'
        %
        controlPolicy
        
    end
    
    
    methods

        function model = ControlModel(inputparams)

            model = model@BaseModel();
            
            fdnames = {'controlPolicy', ...
                       'CRate'};
            model = dispatchParams(model, inputparams, fdnames);
            
        end
        
        function model = registerVarAndPropfuncNames(model)

            model = registerVarAndPropfuncNames@BaseModel(model);
            
            varnames = {};
            
            % Terminal voltage / [V]
            varnames{end + 1} = 'E';
            % Terminal Current / [A]
            varnames{end + 1} = 'I';
            % Equation that relates E and I (depends on geometry and discretization)
            varnames{end + 1} = 'EIequation';
            % control equation
            varnames{end + 1} = 'controlEquation';
            
            model = model.registerVarNames(varnames);

        end

        function cleanState = addStaticVariables(model, cleanState, state)
        % If the control model includes some static variable (that should be initialized at each Newton step), they can be added here.
        % Base class behaviour is do nothing.
        end

        function state = prepareStepControl(model, state, state0, dt, drivingForces)
        % Attach to state the values necessary for the control. This is run only once at the beginning of a time step
        % Base class behaviour is do nothing.
        end
        
        function state = updateControlEquation(model, state)
        % Implemented by child model
        % Base class behaviour is do nothing.
            state.controlEquation = [];
        end
        
        function state = updateControlState(model, state)
        % Implemented by child model.
        % Base class behaviour is do nothing.
        end
        
        function state = updateControlAfterConvergence(model, state, state0, dt)
        % This function is called in updateAfterConvergence after convergence and gives possibility to detect control switch.
        % Base class behaviour is do nothing.
                
        end


        function func = setupStopFunction(model)
        % setup and return a "stop function" for the given control, with signature
        % func = @(model, state, state_prev) ();
        % Default is no stop function (the simulation ends at the last given time step)

            func = @(model, state, state_prev) (false);
            
        end

        function arefulfilled = checkConstraints(model, state)
        % Method that can be implemented to check if the constraints that may be part of the control are fulfilled

        % No default behaviour is implemented here

            error('virtual method');
            
        end


        function func = setupControlFunction(model)
        % Setup and return a "control function". We do not require a special signature for this function. The user
        % should instead make sure the method :code:`updateControl` in  :code:`Battery` cover the implemented control policy.

            error('virtual method');
            
        end

        function step = setupScheduleStep(model, timeSteppingParams)
            
        % Setup and a return the step structure that is part of the schedule which is used as input for
        % :mrst:`simulateScheduleAD`. For some control type, there is a natural construction for this structure. This is
        % why we include this method here, for convenience. It can be overloaded by derived classes. The
        % timeSteppingParams structure by default is given by the data described in :battmofile:`Utilities/JsonSchemas/TimeStepping.schema.json`

            error('virtual method')

            
        end

        function control = setupScheduleControl(model)
        % Setup and return the control structure that is sent to :mrst:`simulateScheduleAD`
            
            control.stopFunction = model.setupStopFunction();
            control.src          = model.setupControlFunction();
            
        end

        function schedule = setupSchedule(model, jsonstruct)
        % Convenience function to setup schedule from main jsonstruct with property TimeStepping
            
            if isfield(jsonstruct, 'TimeStepping')
                timeSteppingParams = jsonstruct.TimeStepping;
            else
                timeSteppingParams = [];
            end

            step    = model.setupScheduleStep(timeSteppingParams);
            control = model.setupScheduleControl();

            schedule = struct('step', step, ...
                              'control', control);
            
        end
        
    end
    
    methods(Static)

        function params = parseTimeSteppingStruct(params)

            paramstemplate = struct('totalTime'          , []   , ...
                                    'numberOfTimeSteps'  , 100  , ...
                                    'timeStepDuration'   , []   , ...
                                    'useRampup'          , false, ...
                                    'numberOfRampupSteps', 5);

            if (nargin > 0) && ~isempty(params)
                
                params = resolveUnitInputJson(params);

                vals    = struct2cell(params);
                fdnames = fieldnames(params);

                params = horzcat(fdnames, vals)';
                params = params(:);


                params = merge_options(paramstemplate, params{:});
                
            else
                
                params = paramstemplate;
                
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
