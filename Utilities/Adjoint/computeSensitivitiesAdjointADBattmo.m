function [sens, lambdas] = computeSensitivitiesAdjointADBattmo(setup, states, params, getObjective, varargin)
% Compute parameter sensitivities using adjoint simulation
%
% SYNOPSIS:
%   sens = computeSensitivitiesAdjointADBattmo(state0, states, model, schedule, getObjective, varargin)
%
% DESCRIPTION:
%   For a given schedule, compute senistivities with regards to parameters
%
% REQUIRED PARAMETERS:
%
%   SimulatorSetup - structure containing:
%
%       state0   - Physical model state at `t = 0`
%       model    - Subclass of PhysicalModel class such as `Battery` that models the physical
%                  effects we want to study.
%       schedule - Schedule suitable for `simulateScheduleAD`.
%
%   states       - All previous states. Must support the syntax
%                  `state = states{i}`. If the problem is too large to fit in
%                  memory, it can use `ResultHandler` class to retrieve files
%                  from the disk.
%
%   params       - cell array of parameters of class ModelParameter
%
%   getObjective - Function handle for getting objective function value
%                  for a given timestep with derivatives. signature: @(tstep, model, state, computeStatePartial)
%
% OPTIONAL PARAMETERS:
%
%   'LinearSolver'   - Subclass of `LinearSolverAD` suitable for solving the
%                      adjoint systems.
%
% RETURNS:
%   sens - Structure with parameter sensitivites of the form
%          sens.(paramName) = paramValue
%

    assert(isa(params{1}, 'ModelParameter'), ...
           'Parameters must be initialized using ''ModelParameter''.')

    opt = struct('LinearSolver', []);
    opt = merge_options(opt, varargin{:});
    if mrstVerbose && setup.model.nonlinearTolerance >= 1e-3
        fprintf(['The accuracy in the gradient depend on the',...
                 'acuracy on the CNV tolerance.\n',...
                 'For better accuracy set a lower value for '...
                 'model.toleranceCNV.'] )
    end

    if isempty(opt.LinearSolver)
        linsolve = BackslashSolverAD();
    else
        linsolve = opt.LinearSolver;
    end

    sens = struct;
    for k = 1 : numel(params)
        sens.(params{k}.name) = 0;
    end
    pNames = fieldnames(sens);

    % Split parameters in inital state/non-initial state due to different  handling
    isInitParam = cellfun(@(p)strcmp(p.belongsTo, 'state0'), params);
    if any(isInitParam)
        error('Modification of the initial state is not supported yet.')
        [initparam, params] = deal(params(isInitParam), params(~isInitParam));
    end

    % inititialize parameters to ADI
    setupParam = initModelParametersADI(setup, params);

    % Propagate AD in model by computing again the parameters in the model that depend on the AD-parameters that have been
    % set directly above.
    modelParam    = setupParam.model;
    scheduleParam = setupParam.schedule;

    modelParam = validateModel(modelParam);

    nstep    = numel(setup.schedule.step.val);
    lambdas   = cell(nstep, 1);
    getState = @(i) getStateFromInput(setup.schedule, states, setup.state0, i);

    getObjectiveState = @(tstep, model, state) getObjective(tstep, model, state, true);
    getObjectiveModel = @(tstep, model, state) getObjective(tstep, model, state, false);

    % Run adjoint
    lambdaVec = [];
    %adjointState = struct();

    for step = nstep : -1 : 1

        fprintf('Solving reverse mode step %d of %d\n', nstep - step + 1, nstep);

        % Compute Lagrange multipliers for the adjoint formulation
        [lambda, lambdaVec]= setup.model.solveAdjoint(linsolve, getState, getObjectiveState, setup.schedule, lambdaVec, step);

        if nargout > 1
            %adjointState = updateAdjointState(setup.model, adjointState, lambda);
            lamdas{step} = lambda;
        end

        % Compute derivatives of the residual equations with respect to the parameters
        [eqdth, objth] = partialWRTparam(modelParam, getState, scheduleParam, step, getObjectiveModel, params);

        % Assemble the sensitivities using the lagrange multipliers
        for kp = 1 : numel(params)
            nm = params{kp}.name;
            for nl = 1 : numel(lambda)
                if isa(eqdth{nl}, 'ADI')
                    sens.(nm) = sens.(nm) + eqdth{nl}.jac{kp}'*lambda{nl};
                end
            end
            if isa(objth, 'ADI')
                sens.(nm) = sens.(nm) + objth.jac{kp}';
            end
        end

    end

    % Compute partial derivative of first eq wrt init state and compute initial state sensitivities
    if any(isInitParam)

        error('Modification of the initial state is not supported yet.')
        schedule = setup.schedule;
        forces   = setup.model.getDrivingForces(schedule.control(schedule.step.control(1)));
        forces   = merge_options(setup.model.getValidDrivingForces(), forces{:});
        model    = setup.model.validateModel(forces);

        state0 = model.validateState(setup.state0);
        % set wellSols just to make subsequent function-calls happy, sensitivities wrt wellSols doesn't make sense anyway
        state0.wellSol = states{1}.wellSol;
        dt = schedule.step.val(1);

        linProblem = model.getAdjointEquations(state0, states{1}, dt, forces,...
                                               'iteration', inf,  ...
                                               'reverseMode', true);

        nms    = applyFunction(@lower, pNames(isInitParam));
        varNms = applyFunction(@lower, linProblem.primaryVariables);

        for k = 1 : numel(nms)
            kn = find(strcmp(nms{k}, varNms));
            assert(numel(kn)==1, 'Unable to match initial state parameter name %s\n', nms{k});
            for nl = 1 : numel(lambda)
                if isa(linProblem.equations{nl}, 'ADI')
                    sens.(nms{k}) = sens.(nms{k}) - linProblem.equations{nl}.jac{kn}'*lambda{nl};
                end
            end
            if strcmp(initparam{k}.type, 'multiplier')
                sens.(nms{k}) = sens.(nms{k}).*initparam{k}.referenceValue;
            end
            sens.(nms{k}) =  initparam{k}.collapseGradient(sens.(nms{k}));
        end

    end

end

function [eqdth, obj] = partialWRTparam(model, getState, schedule, step, getObjective, params)

    validforces = model.getValidDrivingForces();
    current     = getState(step);
    before      = getState(step - 1);
    dt_steps    = schedule.step.val;
    dt          = dt_steps(step);
    cNo         = schedule.step.control(step);

    control = schedule.control(cNo);
    forces  = model.getDrivingForces(control);
    forces  = merge_options(validforces, forces{:});
    model   = model.validateModel(forces);

    % May be not needed
    if step == 1
        before = model.validateState(before);
    end

    % Assemble the equations with AD initialization done such that we obtain the derivative with respect to the parameters
    problem = model.getEquations(before, current, dt, forces, 'iteration', inf, 'resOnly', true);

    eqdth = problem.equations;

    obj = getObjective(step, model, current);
    obj = obj{1};

end

function state = getStateFromInput(schedule, states, state0, i)
    if i == 0
        state = state0;
    elseif i > numel(schedule.step.val)
        state = [];
    else
        state = states{i};
    end
end

function setupParam = initModelParametersADI(setup, param)

    setupParam = setup;

    v  = applyFunction(@(p)p.getParameter(setup), param);
    % use same backend as problem.model
    if isfield(setup, 'model') && isprop(setup.model, 'AutoDiffBackend')
        [v{:}] = setup.model.AutoDiffBackend.initVariablesAD(v{:});
    else
        [v{:}] = initVariablesADI(v{:});
    end

    for k = 1 : numel(v)
        setupParam = param{k}.setParameter(setupParam, v{k});
    end

end

function adjointState = updateAdjointState(model, adjointState, lambda)

    keyboard;

    p = model.getPrimaryVariableNames();
    n = model.equationVarNames();

    for i = 1 : numel(dx)
        val = model.getProp(state, p{i});
        val = val + dx{i};
        state = model.setProp(state, p{i}, val);
    end

    lambdas{step} = lambda;

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
