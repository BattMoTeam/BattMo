function plotBatteryGrid(model, varargin)
%% Plot battery components

    opt = struct('setstyle', true, ...
                 'fig', []);
    opt = merge_options(opt, varargin{:});

    ne    = 'NegativeElectrode';
    pe    = 'PositiveElectrode';
    cc    = 'CurrentCollector';
    sep   = 'Separator';
    co    = 'Coating';

    colors = crameri('vik', 5);
    if isempty(opt.fig)
        figure
    else
        figure(opt.fig);
    end
    hold on
    legtext = {};

    G = model.grid;

    if G.griddim == 1
        edgeparams = {'linewidth', 2};
        facecolorname = 'color';
    else
        edgeparams = {'edgealpha', 0.5, 'edgecolor', 0*[1, 1, 1]};
        facecolorname = 'facecolor';
    end

    if model.include_current_collectors
        plotGrid(model.(pe).(cc).grid, facecolorname, colors(5,:), edgeparams{:});
        legtext{end+1} = 'positive electrode current collector';
    end

    plotGrid(model.(pe).(co).grid, facecolorname, colors(4,:), edgeparams{:});
    plotGrid(model.(sep).grid,     facecolorname, colors(3,:), edgeparams{:});
    plotGrid(model.(ne).(co).grid, facecolorname, colors(2,:), edgeparams{:});
    legtext = [legtext, {'positive electrode active material', 'separator', 'negative electrode active material'}];

    if model.include_current_collectors
        plotGrid(model.(ne).(cc).grid, facecolorname, colors(1,:), edgeparams{:});
        legtext{end+1} = 'negative electrode current collector';
    end

    if G.griddim == 3
        view(3)
    end

    legend(legtext, 'location', 'southwest');
    axis tight;

    if opt.setstyle
        setFigureStyle('quantity', 'single');
    end

    drawnow();

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
