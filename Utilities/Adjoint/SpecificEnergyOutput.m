function obj = SpecificEnergyOutput(model, states, schedule, cutoffVoltage, varargin)
%
%
% SYNOPSIS:
%   function obj = SpecificEnergy(model, states, schedule, varargin)
%
% DESCRIPTION: Computes the specific energy output, sum_{i} ( E_i*I_i*dt_i ) / mass. Here, the index i denotes the time step.
%
%              This function can only be run using 'IEswitch' control where the lowerCutoffVoltage given there is lower than cutoffVoltage
%
% PARAMETERS:
%   model     - Battery model that is used by the solver
%   states    - Input states
%   schedule  - Schedule used for the simulation
%
% KEYWORD ARGUMENTS:
%
%   tStep            - if set, only the given time steps are handled. Otherwise, the whole schedule is used.
%   ComputePartials  - if true, the derivative of the objective functions are also included, see below
%   checkConsistency - if true, run checkConsistency function (this function also serves as documentation)
%
%   
% RETURNS:
%
%   obj - Objective function cell array. One value per time step : obj{i} = E_i*I_i*dt_i. If the option
%         'ComputePartials' is set, the derivative of the objective with
%         respect to state is returned in a format appropriate for the adjoint computation.
%   

    opt = struct('ComputePartials' , false, ...
                 'tStep'           , []   , ...
                 'state'           , []   , ...
                 'from_states'     , true , ...
                 'checkConsistency', false);
    opt     = merge_options(opt, varargin{:});

    if opt.checkConsistency
        assert(checkConsistency(model, schedule, cutoffVoltage));
    end
    
    dts   = schedule.step.val;

    tSteps = opt.tStep;
    
    if isempty(tSteps) % do all
        numSteps = numel(dts);
        tSteps   = (1:numSteps)';
    else
        numSteps = 1;
        dts      = dts(opt.tStep);
    end

    obj = repmat({[]}, numSteps, 1);

    for step = 1:numSteps
        
        dt = dts(step);
        state = states{tSteps(step)};
        if opt.ComputePartials
            if (opt.from_states) 
                state = model.initStateAD(state);
            else
                state = opt.state;
            end
            E = state.Control.E;
            I = state.Control.I;
        else
            E = state.Control.E;
            I = state.Control.I; 
        end

        if value(E) > cutoffVoltage
            obj{step} = I*E*dt;
        else
            obj{step} = double2ADI(0, I*E*dt);
        end
    end

end


function isok = checkConsistency(model, schedule, cutoffVoltage)
    
    isok = true;

    if ~strcmp(model.Control.controlPolicy, 'IEswitch')
        isok = false
        return
    end

    if model.Control.lowerCutoffVoltage > cutoffVoltage
        isok = false
        return
    end

    if numel(schedule.control) > 1 | ~isfield(schedule.control(1), 'IEswitch')
        isok = false;
        return
    end
    
end


%{
Copyright 2021-2023 SINTEF Industry, Sustainable Energy Technology
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
