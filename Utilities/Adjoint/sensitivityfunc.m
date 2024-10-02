function obj = sensitivityfunc(model,states,schedule, varargin)
%
%
% SYNOPSIS:
%   function obj = sensitivityfunc(model, states, refStates, schedule, varargin)
%
% DESCRIPTION: Computes the least-square difference between a reference state and a given computed state for the voltage and current:
%              sum_{i} ( (E_i - Eref_i)^2*dt_i + (I_i - Iref_i)^2*dt_i
%
%
% PARAMETERS:
%   model     - Battery model that is used by the solver
%   states    - Input states
%   refStates - Reference states
%   schedule  - Schedule used for the simulation
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
                 'relTol'         , 1e-10);
    opt = merge_options(opt, varargin{:});

    dts = schedule.step.val;

    if isempty(opt.tStep) %do all
        time     = cumsum(dts);
        numSteps = numel(dts);
    else
        time     = sum(dts(1:(opt.tStep)));
        numSteps = 1;
        dts      = dts(opt.tStep);
    end

    obj = repmat({[]}, numSteps, 1);

    relTol = opt.relTol;

    for k = 1 : numSteps

        t  = time(k);
        dt = dts(k);

        % Find states for time t (given by schedule)
        state    = findState(t, states, dt, relTol);
   

        if opt.ComputePartials
            if opt.from_states
                state = model.initStateAD(state);
            else
                state = opt.state;
            end
        end

        E = state.Control.E;
        I = state.Control.I;

        obj{k} = (E) * dt

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