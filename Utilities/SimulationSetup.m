classdef SimulationSetup
%
%
% SYNOPSIS:
%   SimulationSetup(model, varargin)
%
% DESCRIPTION: Setup a simulation setup instance. The setup of the simulation can be modified by changing the properties in the
%              class. The simulation is run and the output returned using the method `run` (see below)
%              
%              
%
% PARAMETERS: struct with fields
% - model 
% - initstate
% - schedule 
%
% RETURNS:
%   class instance
%
% EXAMPLE: simsetup = SimulationSetup(model);
%          [states, globVars, reports] = simsetup.run();
%
% SEE ALSO: simulateScheduleAD (underlying function that is called for the simulation), 
%

    properties
        
        model % The model that assembles the governing equations, which means,
              % computes the resisuals of the equations and their jacobians. This
              % must be a subclass of the `BaseModel` base class.  Examples are
              % the `Battery`, `GenericBattery` or `Electrolyser` models

        initstate % Initial state of the system. It should have whatever fields
                  % are associated with the physical model, with reasonable
                  % values. It is the responsibility of the user to ensure that
                  % the state is properly initialized.
                  %
                  % If the model has a `setupInitialState` method, this method
                  % is used to setup the `initialState` (if none is given)
                  % 

        schedule  % Schedule containing fields step and control, defined as follows:
                  % - `schedule.control` is a struct array containing
                  %    fields that the model knows how to process.
                  %    Typically, this will be the fields such as `.src` 
                  %
                  % - `schedule.step` contains two arrays of equal 
                  %    size named `val` and `control`. Control is a
                  %    index into the `schedule.control` array,
                  %    indicating which control is to be used for the
                  %    timestep.`schedule.step.val` is the timestep
                  %    used for that control step.
                  %
                  % If the model has a `Control` submodel which provides a
                  % `setupSchedule` method, then this method is used to setup the
                  % `schedule` (if none is given)
                  % 
        
        Verbose % Indicate if extra output is to be printed such as detailed
                % convergence reports and so on.

        OutputMinisteps % The solver may not use timesteps equal to the control
                        % steps depending on problem stiffness and timestep
                        % selection. Enabling this option will make the solver
                        % output the states and reports for all steps actually
                        % taken and not just at the control steps. See also
                        % `convertReportToSchedule` which can be used to construct
                        % a new schedule from these timesteps. Default value is true.

        NonLinearSolver  % An instance of the `NonLinearSolver` class. Consider
                         % using this if you for example want a special timestep
                         % selection algorithm. See the `NonLinearSolver` class
                         % docs for more information.

        OutputHandler    % Output handler class, for example for writing states
                         % to disk during the simulation or in-situ
                         % visualization. See the ResultHandler base class.

        GlobVarOutputHandler % Same as OutputHandler, but for the well
                              % solutions for the individual report steps. Well
                              % solutions are also stored using OutputHandler, but
                              % using GlobVarOutputHandler is convenient for quickly
                              % loading well solutions only.

        ReportHandler % Same as OutputHandler, but for the reports for the report steps.

        LinearSolver % Class instance subclassed from `LinearSolverAD`. Used to
                     % solve linearized problems in the `NonLinearSolver`
                     % class. Note that if you are passing a `NonLinearSolver`,
                     % you might as well put it there.

        afterStepFn    % Function handle to an optional function that will be
                       % called after each successful report step in the
                       % schedule. The function should take in the following
                       % input arguments:
                       %  - model: The model used in the schedule
                       %  - states: A cell array of all states that are
                       %    computed, as well as possible empty entries
                       %    where the states have not been computed yet.
                       %  - reports: A cell array of reports for each step,
                       %    with empty entries for steps that have not been
                       %    reached yet.
                       %  - solver: The NonLinearSolver instance.
                       %  - schedule: The current schedule.
                       %  - simtime: Array with the time in seconds taken by
                       %    the `NonLinearSolver` to compute each step.
                       %    Entries not computed will contain zeros.
                       %
                       %    See `getPlotAfterStep` for more information and
                       %    `blackoilTutorialPlotHook` for a worked example.

        processOutputFn  % Function handle to an optional function that
                         % processes the simulation output (globVars, states
                         % and reports) before they are stored to state using
                         % the output handlers. Allows for storing only data
                         % of interest, which is useful when dealing with
                         % large models. Changes made to the output by this
                         % function are only applied to the data that is
                         % stored, and will not affect what is passed on to
                         % the next timestep.
   
        controlLogicFn  % - Function handle to optional function that will be
                        %     called after each step enabling schedule updates to 
                        %     be triggered on specified events. Input arguemnts:
                        %     - state: The current state
                        %     - schedule: The current schedule
                        %     - report: Current report
                        %     - i: The current report step such that current
                        %       control step equals schedule.step.control(i)
                        %     The function must have three outputs:
                        %     - schedule: Possibly updated schedule
                        %     - report: Possibly updated report
                        %     - isAltered: Flag indicated whether the schedule
                        %       was updated or not

    end

    methods
        
        function simsetup = SimulationSetup(simInput)

            model = simInput.model;
            
            simsetup.model = model;

            %% Setup default initial state, if possible

            if isUnAssigned(simInput, 'initstate')
                if ismethod(model, 'setupInitialState')
                    simInput = setJsonStructField(simInput, 'initstate', model.setupInitialState());
                end
            end
            
            %% Setup default schedule (time stepping), if possible
            
            if isUnAssigned(simInput, 'schedule')
                if isAssigned(model, 'Control') && ismethod(model.Control, 'setupSchedule')
                    simInput = setJsonStructField(simInput, 'schedule', model.Control.setupSchedule());
                end
            end
            
            %% Set default nonlinear solver
            
            nls = NonLinearSolver();
            nls.errorOnFailure    = false;
            nls.continueOnFailure = false;

            simInput = setDefaultJsonStructField(simInput, 'NonLinearSolver', nls);

            simsetup.schedule        = getJsonStructField(simInput, 'schedule');
            simsetup.initstate       = getJsonStructField(simInput, 'initstate');
            simsetup.NonLinearSolver = getJsonStructField(simInput, 'NonLinearSolver');
            simsetup.OutputMinisteps = getJsonStructField(simInput, 'OutputMinisteps', true);
            
        end

        function [states, globVars, reports] = run(simsetup)

            initstate = simsetup.initstate;
            model     = simsetup.model;
            schedule  = simsetup.schedule;

            fds = {'Verbose'             , ...
                   'OutputMinisteps'     , ...
                   'Verbose'             , ...
                   'OutputMinisteps'     , ...
                   'NonLinearSolver'     , ...
                   'OutputHandler'       , ...
                   'GlobVarOutputHandler', ...
                   'ReportHandler'       , ...
                   'LinearSolver'        , ...
                   'afterStepFn'         , ...
                   'processOutputFn'     , ...
                   'controlLogicFn'};

            vals = cellfun(@(fd) simsetup.(fd), fds, 'uniformoutput', false);

            opts = reshape(vertcat(fds, vals), [], 1);
            
            [globVars, states, reports] = simulateScheduleAD(initstate, model, schedule, opts{:});
            
        end
        
    end
    
end



