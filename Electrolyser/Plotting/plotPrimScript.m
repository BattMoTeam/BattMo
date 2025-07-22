setupoutput = false;
if setupoutput
    ind = cellfun(@(state) ~isempty(state), states);
    states = states(ind);
    for istate = 1 : numel(states)
        states{istate} = model.addVariables(states{istate});
    end
    time = cellfun(@(state) state.time, states);
    E = cellfun(@(state) state.(oer).(ptl).E, states);
    I = cellfun(@(state) state.(oer).(ctl).I, states);

    close all

    figure
    plot(time/hour, E)
    xlabel('time [hour]');
    ylabel('voltage');

    figure
    plot(time/hour, -I)
    xlabel('time [hour]');
    ylabel('Current [A]');
end

close all

fds = {'HydrogenEvolutionElectrode.PorousTransportLayer.H2rhoeps', ...
       'HydrogenEvolutionElectrode.PorousTransportLayer.liqrhoeps', ...
       'HydrogenEvolutionElectrode.PorousTransportLayer.liqeps', ...
       'HydrogenEvolutionElectrode.PorousTransportLayer.H2Ogasrhoeps', ...
       'HydrogenEvolutionElectrode.PorousTransportLayer.OHceps', ...
       'HydrogenEvolutionElectrode.PorousTransportLayer.phi', ...
       'OxygenEvolutionElectrode.PorousTransportLayer.O2rhoeps', ...
       'OxygenEvolutionElectrode.PorousTransportLayer.liqrhoeps', ...
       'OxygenEvolutionElectrode.PorousTransportLayer.liqeps', ...
       'OxygenEvolutionElectrode.PorousTransportLayer.H2Ogasrhoeps', ...
       'OxygenEvolutionElectrode.PorousTransportLayer.OHceps', ...
       'OxygenEvolutionElectrode.PorousTransportLayer.phi', ...
       'IonomerMembrane.H2Oceps', ...
       'IonomerMembrane.phi'};

dosave = true;

nt = numel(time);

for ivar = 1 : numel(fds)
    fd = fds{ivar};

    varind = cgti.regexpVarNameSelect(fd);
    assert(numel(varind) == 1, ['too many field selected for ' fd]);
    varname = cgti.varNameList{varind};
    nodenames = cgti.getNodeName(varname);
    % should be unique here
    nodename = nodenames{1};
    [ymin, ymax, y] = getvals(nodename, states);

    figure
    hold on
    dominmax = false;
    if dominmax
        plot(time, ymin);
        plot(time, ymax);
    else
        imax = numel(time);
        for itime = 1 : imax
            h = plot(y{itime});
            if itime == 1
                set(h, 'linewidth', 3, 'color', 'b');
            elseif itime == imax
                set(h, 'linewidth', 3, 'color', 'r');
                title(nodename, ' ');
                if dosave
                    filename = sprintf('%s.png', nodename);
                    filename = fullfile(battmoDir, 'Electrolyser','img', filename);
                    saveas(gcf, filename);
                end
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
