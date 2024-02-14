function  [model, states, reports, solver, ok] = plotAfterStepIV(model, states, reports, solver, schedule, simtime)
%% Process output and recover the output voltage and current from the output states.

    ind = cellfun(@(x) not(isempty(x)), states);
    states = states(ind);

    E    = cellfun(@(x) x.Control.E, states);
    I    = cellfun(@(x) x.Control.I, states);
    time = cellfun(@(x) x.time, states);

    dt = diff([0, time]);
    SOC = cumsum(I.*dt);

    %%
    if ~ishandle(33)
        figure(33)
    else
        % Avoid stealing focus if figure already exists
        set(0, 'CurrentFigure', 33);
    end

    subplot(2, 2, 1)
    hold on; grid on
    plot(time/hour, I, '*-')
    xlabel('time  / h');
    ylabel('current  / A');

    subplot(2, 2, 2),
    hold on; grid on
    plot(SOC/hour, E, '*-')
    xlabel('SOC  / Ah');
    ylabel('voltage  / V');

    subplot(2, 2, 3),
    hold on; grid on
    plot(time/hour, E.*I.*dt/hour, '*-')
    xlabel('time  / h');
    ylabel('power  / W');

    subplot(2, 2, 4),
    hold on; grid on
    plot(SOC, cumsum(E.*I.*dt)/hour, '*-')
    xlabel('SOC  / Ah');
    ylabel('energy  / Wh');

    %% Temperature

    if model.use_thermal
        Tmax = cellfun(@(x) max(x.ThermalModel.T), states);

        figure(44)
        hold on; grid on

        subplot(2, 1, 1)
        plot(time/hour, Tmax, '*-')
        xlabel('time  / h');
        ylabel('maximum temperature  / K');

        subplot(2, 1, 2)
        plot(time/hour, cumsum(I.*dt)/hour, '*-')
        xlabel('time  / h');
        ylabel('capacity  / Ah');
    end

    ok = true;
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
