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

% load results from old code

close all

% str = 'H P ressure';
% str = ' P conc 2';
% str = 'H P j$';
str = 'P phi$';
% str = 'Cat Lay phi';
% str = 'phi$';
% str = 'cOHelyte';
% str = 'Cat Rate';
% str = 'Ox P B cOH';

% cgit.printVarNames(str);

% varinds = cgit.regexpVarNameSelect(str);

% fds = {'OxygenEvolutionElectrode.PorousTransportLayer.phasePressure 1', ...
%        'OxygenEvolutionElectrode.PorousTransportLayer.phasePressure 2', ...
%        'HydrogenEvolutionElectrode.PorousTransportLayer.phasePressure 1', ...
%        'HydrogenEvolutionElectrode.PorousTransportLayer.phasePressure 2'};

fds = {'OxygenEvolutionElectrode.PorousTransportLayer.concentration 1'  , ...
       'OxygenEvolutionElectrode.PorousTransportLayer.concentration 2'  , ...
       'OxygenEvolutionElectrode.PorousTransportLayer.concentration 3'  , ...
       'HydrogenEvolutionElectrode.PorousTransportLayer.concentration 1', ...
       'HydrogenEvolutionElectrode.PorousTransportLayer.concentration 2', ...
       'HydrogenEvolutionElectrode.PorousTransportLayer.concentration 3'};

% fds = {'OxygenEvolutionElectrode.PorousTransportLayer.H2OvaporLiquidExchangeRate'  , ...
%        'OxygenEvolutionElectrode.ExchangeReaction.H2OexchangeRate'                    , ...
%        'OxygenEvolutionElectrode.ExchangeReaction.OHexchangeRate'                     , ...
%        'HydrogenEvolutionElectrode.PorousTransportLayer.H2OvaporLiquidExchangeRate', ...
%        'HydrogenEvolutionElectrode.ExchangeReaction.H2OexchangeRate'                  , ...
%        'HydrogenEvolutionElectrode.ExchangeReaction.OHexchangeRate'};

% fds = {'HydrogenEvolutionElectrode.PorousTransportLayer.H2Oa', ...
%        'OxygenEvolutionElectrode.PorousTransportLayer.H2Oa', ...
%        'IonomerMembrane.H2Oa$'};

% fds = {'OxygenEvolutionElectrode.PorousTransportLayer.liqrho$', ...
       % 'HydrogenEvolutionElectrode.PorousTransportLayer.liqrho$'};

% fds = {'OxygenEvolutionElectrode.PorousTransportLayer.volumeFr 1', ...
%        'OxygenEvolutionElectrode.PorousTransportLayer.volumeFr 2', ...
%        'HydrogenEvolutionElectrode.PorousTransportLayer.volumeFr 1', ...
%        'HydrogenEvolutionElectrode.PorousTransportLayer.volumeFr 2'};

% fds = {'OxygenEvolutionElectrode.PorousTransportLayer.OHsource', ...
%        'HydrogenEvolutionElectrode.PorousTransportLayer.OHsource'};


dosave = false;

nt = numel(time);

for ivar = 1 : numel(fds)
    fd = fds{ivar};

    varind = cgit.regexpVarNameSelect(fd);
    assert(numel(varind) == 1, ['too many field selected for ' fd]);
    varname = cgit.varNameList{varind};
    nodenames = cgit.getNodeName(varname);
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
                title(nodenames, ' ');
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
