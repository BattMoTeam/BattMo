function obj = leastSquaresEI(simsetup, states, refStates, varargin)
%
%
% SYNOPSIS:
%   function obj = leastSquaresEI(model, states, refStates, schedule, varargin)
%
% DESCRIPTION: Computes the least-square difference between a reference state and a given computed state for the voltage and current:
%              sum_{i} ( (E_i - Eref_i)^2*dt_i + (I_i - Iref_i)^2*dt_i
%
%
% PARAMETERS:
%   simsetup  - Instance of SimulatorSetup
%   states    - Input states
%   refStates - Reference states
%
% KEYWORD ARGUMENTS:
%
%   tStep           - if set, only the given time steps are handled. Otherwise, the whole schedule is used.
%   ComputePartials - if true, the derivative of the objective functions are also included, see below
%
%
% RETURNS:
%
%   obj   - Objective function cell array. One value per time step : obj{i} = (E_i - Eref_i)^2*dt_i + (I_i -
%           Iref_i)^2*dt_i. If the option 'ComputePartials' is set, the derivative of the objective with
%           respect to state is returned in a format appropriate for the adjoint computation.
%

    opt = struct('ComputePartials', false, ...
                 'tStep'          , []   , ...
                 'state'          , []   , ...
                 'from_states'    , true , ...
                 'relTol'         , 1e-10, ...
                 'scaling'        , []   , ...
                 'findState'      , false, ...
                 'includeI'       , false);

    opt = merge_options(opt, varargin{:});

    model    = simsetup.model;
    schedule = simsetup.schedule;
    
    dts = schedule.step.val;

    if isempty(opt.tStep) % do all
        numSteps = numel(dts);
        time     = cumsum(dts);
        tSteps   = (1:numSteps)';
    else
        numSteps = 1;
        dts      = dts(opt.tStep);
        tSteps   = opt.tStep;
    end

    obj = repmat({[]}, numSteps, 1);

    relTol = opt.relTol;

    for k = 1 : numSteps

        dt = dts(k);
        tStep = tSteps(k);

        if opt.findState
            % Find states for time t (given by schedule)
            t        = time(k);
            state    = findState(t, states, dt, relTol);
            stateRef = findState(t, refStates, dt, relTol);
        else
            state    = states{tStep};
            stateRef = refStates{tStep};
        end

        if opt.ComputePartials
            if opt.from_states
                state = model.initStateAD(state);
            else
                state = opt.state;
            end
        end

        E      = state.Control.E;
        Eref   = stateRef.Control.E;
        obj{k} = (E - Eref)^2 * dt;

        if opt.includeI
            I      = state.Control.I;
            Iref   = stateRef.Control.I;
            obj{k} = obj{k} + (I - Iref)^2 * dt;
        end

    end

    if ~isempty(opt.scaling)
        obj = applyFunction(@(x) x/opt.scaling, obj);
    end

end


function state = findState(tstar, states, dt, relTol)

% Find state with time matching tstar

    times = cellfun(@(x) x.time, states);
    [tdiff, idx] = min(abs(times - tstar));

    assert(tdiff / dt < relTol);

    state = states{idx};

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
