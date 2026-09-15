function output = setupSimulationForJuliaBridge(jsonstruct, varargin)

    opt = struct('runSimulation', false);
    opt = merge_options(opt, varargin{:});
    
    output = runBattery(jsonstruct, 'runSimulation', opt.runSimulation);

    model     = output.model;
    schedule  = output.schedule;

    if opt.runSimulation

        output.use_state_ref = true;
        
        states   = output.states;

        %% added all the extra variables on state

        for istate = 1 : numel(states)
            states{istate} = model.addVariables(states{istate});
        end

        output.states = states;
    else

        output.use_state_ref = false;        
        
    end

    model = convertModelForJuliaBridge(model);

    model    = class2data(model);
    schedule = class2data(schedule);

    output.model    = model;
    output.schedule = schedule;
    
end

%{
Copyright 2021-2026 SINTEF Industry, Sustainable Energy Technology
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
