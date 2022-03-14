function [] = plotChemical(model, states, varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


%% Parse inputs
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

%% Instantiate figure(s)
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
        % Plot Charge Carrier Source
        subplot(1,2,1), plotChargeCarrierSource(model,states{i})
            colormap(gca, fs.colormap_concentration);
            %caxis([min(states{1}.Electrolyte.LiSource), max(states{end}.Electrolyte.LiSource)]);
            colorbar();
            cbar = get(gca, 'Colorbar');
            cbar.Color = fs.fontColor;
            title('Heat Source  /  mol \cdot m^{-3} \cdot s^{-1}', 'color', fs.fontColor)
            set(gca, 'XColor', fs.fontColor);
            set(gca, 'YColor', fs.fontColor);
            
        % Plot Charge Carrier Concentration
        subplot(1,2,2), plotChargeCarrierConcentration(model,states{i})
            colormap(gca, fs.colormap_concentration);
            %caxis([min(states{1}.Electrolyte.cs{1}), max(states{end}.Electrolyte.cs{1})]);
            colorbar();
            cbar = get(gca, 'Colorbar');
            cbar.Color = fs.fontColor;
            title('Concentration  /  mol \cdot m^{-3}', 'color', fs.fontColor)
            set(gca, 'XColor', fs.fontColor);
            set(gca, 'YColor', fs.fontColor);
        else
            
            % Plot Charge Carrier Source
            figure(1), plotChargeCarrierSource(model,states{i})
            colormap(gca, fs.colormap_concentration);
            caxis([min(states{1}.Electroltye.LiSource), max(states{end}.Electroltye.LiSource)]);
            colorbar();
            cbar = get(gca, 'Colorbar');
            cbar.Color = fs.fontColor;
            title('Heat Source  /  mol \cdot m^{-3} \cdot s^{-1}', 'color', fs.fontColor)
            set(gca, 'XColor', fs.fontColor);
            set(gca, 'YColor', fs.fontColor);
            
            % Plot Charge Carrier Concentration
            figure(2), plotChargeCarrierConcentration(model,states{i})
            colormap(gca, fs.colormap_concentration);
            caxis([min(states{1}.Electrolyte.cs{1}), max(states{end}.Electrolyte.cs{1})]);
            colorbar();
            cbar = get(gca, 'Colorbar');
            cbar.Color = fs.fontColor;
            title('Concentration  /  mol \cdot m^{-3}', 'color', fs.fontColor)
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
