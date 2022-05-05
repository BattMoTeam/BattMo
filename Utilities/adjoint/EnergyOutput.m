function obj = EnergyOutput(model, states, schedule, varargin)
% Compute net present value of a schedule with well solutions

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
opt     = struct('Price',             1.0 , ...
                 'ComputePartials',      false, ...
                 'tStep' ,               [],   ...
                 'state',                 [], ...
                 'from_states',          true);
opt     = merge_options(opt, varargin{:});

% pressure and saturaton vectors just used for place-holding
%p  = zeros(G.cells.num, 1);
%sW = zeros(G.cells.num, 1);

dts   = schedule.step.val;

tSteps = opt.tStep;
if isempty(tSteps) %do all
    time = 0;
    numSteps = numel(dts);
    tSteps = (1:numSteps)';
else
    time = sum(dts(1:(opt.tStep-1)));
    numSteps = 1;
    dts = dts(opt.tStep);
end

obj = repmat({[]}, numSteps, 1);

for step = 1:numSteps
    dt = dts(step);
    state = states{tSteps(step)};
    if opt.ComputePartials
         %[~, ~, qWs, qOs, bhp] = ...
        %  initVariablesADI(p, sW, qWs, qOs, bhp);
        if(opt.from_states) 
            %init=true;
            %state = model.getStateAD( states{tSteps(step)}, init);
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
    obj{step} = opt.Price*I*E*dt;
end
end