function [fig] = plotElectrolyteChargeCarrierConcentration(model, states, varargin)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

%
close all

%% Parse inputs            
defaultStep         = length(states);
defaultTheme        = 'light';
expectedTheme       = {'dark', 'light', 'blue', 'offwhite'};
defaultSize         = 'A4';
expectedSize        = {'A4', 'square', 'wide'};
defaultOrientation  = 'landscape';
expectedOrientation = {'landscape', 'portrait'};

p = inputParser;
validStep = @(x) isnumeric(x) && isinteger(int8(x)) && (x >= 0) && (x <= length(states));
addParameter(p, 'step', defaultStep, validStep);
addParameter(p, 'theme', defaultTheme, ...
    @(x) any(validatestring(x,expectedTheme)));
addParameter(p, 'size', defaultSize, ...
    @(x) any(validatestring(x,expectedSize)));
addParameter(p, 'orientation', defaultOrientation, ...
    @(x) any(validatestring(x,expectedOrientation)));
parse(p, varargin{:});

step = p.Results.step;

fig = figure;

if step ~= 0
    if model.G.griddim == 1
        plotCellData(model.Electrolyte.G, states{step}.Electrolyte.c ./ 1000);
        xlabel(gca, 'Position  /  m')
        ylabel(gca, 'Concentration  /  mol \cdot L^{-1}')

    elseif model.G.griddim == 2
        plotCellData(model.Electrolyte.G, states{step}.Electrolyte.c ./ 1000, 'edgealpha', 0.1);
        colormap(cmocean('curl'));
        xlabel(gca, 'Position  /  m')
        ylabel(gca, 'Position  /  m')
        title('Concentration  /  mol \cdot L^{-1}')
        cmin = min(states{step}.Electrolyte.c ./ 1000);
        cmax = max(states{step}.Electrolyte.c ./ 1000);
        cnom = mean(states{1}.Electrolyte.c ./ 1000);
        deltac = max(abs(cmin-cnom), abs(cmax-cnom));
        cmin = cnom-deltac;
        cmax = cnom+deltac;
        caxis([cmin, cmax]);
        colorbar;

    elseif model.G.griddim == 3
        plotCellData(model.Electrolyte.G, states{step}.Electrolyte.c ./ 1000, 'edgealpha', 0.1);
        colormap(cmocean('curl'));
        xlabel(gca, 'Position  /  m')
        ylabel(gca, 'Position  /  m')
        zlabel(gca, 'Position  /  m')
        title('Concentration  /  mol \cdot L^{-1}')
        cmin = min(states{step}.Electrolyte.c ./ 1000);
        cmax = max(states{step}.Electrolyte.c ./ 1000);
        cnom = mean(states{1}.Electrolyte.c ./ 1000);
        deltac = max(abs(cmin-cnom), abs(cmax-cnom));
        cmin = cnom-deltac;
        cmax = cnom+deltac;
        caxis([cmin, cmax]);
        colorbar;
        view(45,45)
    end
    setFigureStyle('theme', p.Results.theme, 'size', p.Results.size, 'orientation', p.Results.orientation, 'quantity', 'single');
else
    for i = 1:length(states)
        if i == 1
            cmax = max(max(states{i}.Electrolyte.c ./ 1000));
            cmin = min(max(states{i}.Electrolyte.c ./ 1000));
            xmin = min(model.Electrolyte.G.nodes.coords(:,1));
            xmax = max(model.Electrolyte.G.nodes.coords(:,1));
            if model.G.griddim == 2
                ymin = min(model.Electrolyte.G.nodes.coords(:,2));
                ymax = max(model.Electrolyte.G.nodes.coords(:,2));
            elseif model.G.griddim == 3
                ymin = min(model.Electrolyte.G.nodes.coords(:,2));
                ymax = max(model.Electrolyte.G.nodes.coords(:,2));
                zmin = min(model.Electrolyte.G.nodes.coords(:,3));
                zmax = max(model.Electrolyte.G.nodes.coords(:,3));
            end
        else
            cmax = max( cmax, max(max(states{i}.Electrolyte.c ./ 1000)));
            cmin = min( cmin, min(min(states{i}.Electrolyte.c ./ 1000)));
        end
    end
    
    for i = 1:length(states)
        if i == 1
            style = setFigureStyle('theme', p.Results.theme, 'size', p.Results.size, 'orientation', p.Results.orientation, 'quantity', 'single'); 
        end
        if model.G.griddim == 1
            plotCellData(model.Electrolyte.G, states{i}.Electrolyte.c ./ 1000, 'LineWidth', style.lineWidth);
            xlabel(gca, 'Position  /  m', 'FontSize', style.fontSize)
            ylabel(gca, 'Concentration  /  mol \cdot L^{-1}', 'FontSize', style.fontSize)
            set(gca, ...
                'FontSize', style.fontSize, ...
                'FontName', style.fontName, ...
                'color', style.backgroundColor, ...
                'XColor', style.fontColor, ...
                'YColor', style.fontColor, ...
                'GridColor', style.fontColor)
            ylim([cmin, cmax]);
            xlim([xmin, xmax]);

        elseif model.G.griddim == 2
            plotCellData(model.Electrolyte.G, states{i}.Electrolyte.c ./ 1000, 'edgealpha', 0.1);
            colormap(cmocean('curl'));
            xlabel(gca, 'Position  /  m')
            ylabel(gca, 'Position  /  m')
            title('Concentration  /  mol \cdot L^{-1}')
            
            xlim([xmin, xmax]);
            ylim([ymin, ymax]);
            cnom = mean(states{1}.Electrolyte.c ./ 1000);
            deltac = max(abs(cmin-cnom), abs(cmax-cnom));
            cmin = cnom-deltac;
            cmax = cnom+deltac;
            caxis([cmin, cmax]);
            c = colorbar;
            c.Color = style.fontColor;

        elseif model.G.griddim == 3
            plotCellData(model.Electrolyte.G, states{i}.Electrolyte.c ./ 1000, 'edgealpha', 0.1);
            colormap(cmocean('curl'));
            xlabel(gca, 'Position  /  m')
            ylabel(gca, 'Position  /  m')
            zlabel(gca, 'Position  /  m')
            title('Concentration  /  mol \cdot L^{-1}')
            xlim([xmin, xmax]);
            ylim([ymin, ymax]);
            zlim([zmin, zmax]);
            cnom = mean(states{1}.Electrolyte.c ./ 1000);
            deltac = max(abs(cmin-cnom), abs(cmax-cnom));
            cmin = cnom-deltac;
            cmax = cnom+deltac;
            caxis([cmin, cmax]);
            c = colorbar;
            c.Color = style.fontColor;
            view(45,45)
        end
        
        drawnow
        pause(0.1)
    end
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
