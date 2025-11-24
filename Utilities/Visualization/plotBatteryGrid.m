function fig = plotBatteryGrid(model, varargin)
%% Plot battery components

    opt = struct('setstyle', true, ...
                 'fig', [], ...
                 'shortLegendText', false, ...
                 'legendLocation', 'sw', ...
                 'axisLabels', false, ...
                 'legend', true);
    opt = merge_options(opt, varargin{:});

    ne    = 'NegativeElectrode';
    pe    = 'PositiveElectrode';
    cc    = 'CurrentCollector';
    sep   = 'Separator';
    co    = 'Coating';

    colors = crameri('vik', 5);
    if isempty(opt.fig)
        fig = figure;
    else
        fig = figure(opt.fig);
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
        legtext{end+1} = 'Positive electrode current collector';
    end

    plotGrid(model.(pe).(co).grid, facecolorname, colors(4,:), edgeparams{:});
    plotGrid(model.(sep).grid,     facecolorname, colors(3,:), edgeparams{:});
    plotGrid(model.(ne).(co).grid, facecolorname, colors(2,:), edgeparams{:});
    legtext = [legtext, {'Positive electrode active material', 'Separator', 'Negative electrode active material'}];

    if model.include_current_collectors
        plotGrid(model.(ne).(cc).grid, facecolorname, colors(1,:), edgeparams{:});
        legtext{end+1} = 'Negative electrode current collector';
    end

    if G.griddim == 3
        view(3)
    end

    if opt.axisLabels
        xlabel 'x';
        if G.griddim > 1
            ylabel 'y';
            if G.griddim > 2
                zlabel 'z';
            end
        end
    end

    if opt.shortLegendText
        legtext = cellfun(@(s) strrep(s, 'Positive electrode', 'PE'), legtext, 'un', false);
        legtext = cellfun(@(s) strrep(s, 'Negative electrode', 'NE'), legtext, 'un', false);
    end

    if opt.legend
        legend(legtext, 'location', opt.legendLocation);
    end
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
