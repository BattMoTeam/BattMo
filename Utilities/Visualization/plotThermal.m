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
addOptional(p, 'docked', defaultDocked, ...
    @(x) any(validatestring(x,expectedDocked)));
addOptional(p, 'theme', defaultTheme, ...
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

