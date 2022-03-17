function [] = plotThermal(model, states, varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

defaultDocked       = 'true';
expectedDocked      = {'true', 'false'};
defaultTheme        = 'dark';
expectedTheme       = {'dark', 'light', 'blue'};

p = inputParser;
addRequired(p, 'model')
addRequired(p, 'states')
addParameter(p, 'docked', defaultDocked, ...
    @(x) any(validatestring(x,expectedDocked)));
addParameter(p, 'theme', defaultTheme, ...
    @(x) any(validatestring(x,expectedTheme)));
parse(p, model, states, varargin{:});

%dimensionality = model.G.griddim;

if strcmpi(p.Results.docked, 'true')
    nFigs = 1;
else
    nFigs = 2;
end

% Set figure style
fs = FigureStyle('theme', p.Results.theme);
figureHandle = cell(nFigs);
for n = 1:nFigs
    figureHandle{n} = figure(n);
    fs.applyFigureStyle(figureHandle{n});
end

    for i = 1:length(states)
        if strcmpi(p.Results.docked, 'true')
        % Plot heat source
        subplot(1,2,1), plotHeatSource(model,states{i})
            colormap(gca, fs.colormap_heatSource);
            caxis([min(states{1}.ThermalModel.jHeatSource), max(states{end}.ThermalModel.jHeatSource)]);
            colorbar();
            cbar = get(gca, 'Colorbar');
            cbar.Color = fs.fontColor;
            title('Heat Source  /  I \cdot m^{-3}', 'color', fs.fontColor)
            set(gca, 'XColor', fs.fontColor);
            set(gca, 'YColor', fs.fontColor);
            
        % Plot Temperature
        subplot(1,2,2), plotTemperature(model,states{i})
            colormap(gca, fs.colormap_temperature);
            caxis([min(states{1}.ThermalModel.T), max(states{end}.ThermalModel.T)]);
            colorbar();
            cbar = get(gca, 'Colorbar');
            cbar.Color = fs.fontColor;
            title('Temperature  /  K', 'color', fs.fontColor)
            set(gca, 'XColor', fs.fontColor);
            set(gca, 'YColor', fs.fontColor);
        else
            
            % Plot Heat Source
            figure(1), plotHeatSource(model,states{i})
            colormap(gca, fs.colormap_heatSource);
            caxis([min(states{1}.ThermalModel.jHeatSource), max(states{end}.ThermalModel.jHeatSource)]);
            colorbar();
            cbar = get(gca, 'Colorbar');
            cbar.Color = fs.fontColor;
            title('Heat Source  /  I \cdot m^{-3}', 'color', fs.fontColor)
            set(gca, 'XColor', fs.fontColor);
            set(gca, 'YColor', fs.fontColor);
            
            % Plot Temperature
            figure(2), plotTemperature(model,states{i})
            colormap(gca, fs.colormap_temperature);
            caxis([min(states{1}.ThermalModel.T), max(states{end}.ThermalModel.T)]);
            colorbar();
            cbar = get(gca, 'Colorbar');
            cbar.Color = fs.fontColor;
            title('Temperature  /  K', 'color', fs.fontColor)
            set(gca, 'XColor', fs.fontColor);
            set(gca, 'YColor', fs.fontColor);
        end
        drawnow
    end

end




%{
Copyright 2021-2022 SINTEF Industry, Sustainable Energy Technology
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
