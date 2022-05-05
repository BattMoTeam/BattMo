function style = setFigureStyle(varargin)
%SETFIGURESTYLE Sets the style of the open figures according to some
%pre-defined templates.
%   Detailed explanation goes here


%% Parse inputs            
defaultTheme        = 'light';
expectedTheme       = {'dark', 'light', 'blue', 'offwhite'};
defaultSize         = 'A4';
expectedSize        = {'A4', 'square', 'wide'};
defaultOrientation  = 'landscape';
expectedOrientation = {'landscape', 'portrait'};
defaultQnty         = 'all';
expectedQnty        = {'single', 'all'};

p = inputParser;
addParameter(p, 'theme', defaultTheme, ...
    @(x) any(validatestring(x,expectedTheme)));
addParameter(p, 'size', defaultSize, ...
    @(x) any(validatestring(x,expectedSize)));
addParameter(p, 'orientation', defaultOrientation, ...
    @(x) any(validatestring(x,expectedOrientation)));
addParameter(p, 'quantity', defaultQnty, ...
    @(x) any(validatestring(x,expectedQnty)));
parse(p, varargin{:});

if strcmpi(p.Results.quantity,'all')
    h = findobj('type','figure');
    n = length(h);
else
    n = 1;
end

for fig = 1 : n

    set(0, 'CurrentFigure', fig);
    ax = gca;
    hLegend = findobj(gcf, 'Type', 'Legend');
    hColorbar = findobj(gcf, 'Type', 'Colorbar');

    %% Define font properties
    style.fontName = 'Helvetica';
    style.fontSize = 24;

    %% Define size properties
    if strcmpi(p.Results.size, 'A4')
        x = 21;    % cm
        x2 = 29.7;  % cm
        if strcmpi(p.Results.orientation, 'landscape')
            style.width   = max(x, x2);
            style.height  = min(x, x2);
        elseif strcmpi(p.Results.orientation, 'portrait')
            style.width   = min(x, x2);
            style.height  = max(x, x2);
        end
    elseif strcmpi(p.Results.size, 'square')
        style.width = 10;    % cm
        style.height = 10;  % cm
    elseif strcmpi(p.Results.size, 'wide')
        style.width = 50.8000;    % cm
        style.height =  25.4794;  % cm
    end

    style.lineWidth = 5;

    %% Define color properties for standard plot themes
    if strcmpi(p.Results.theme,'dark')
        style.backgroundColor         = coal();
        style.fontColor               = 'w';
    elseif strcmpi(p.Results.theme,'light')
        style.backgroundColor         = 'w';
        style.fontColor               = 'k';
    elseif strcmpi(p.Results.theme,'blue')
        style.backgroundColor         = sintefblue();
        style.fontColor               = 'w';
    elseif strcmpi(p.Results.theme,'offwhite')
        style.backgroundColor         = offwhite();
        style.fontColor               = darkgrey();
    end

    %% Apply properties
    set(fig, ...
        'units','centimeter','position',[1, 1, style.width, style.height,], ...
        'color', style.backgroundColor)

    set(ax, ...
        'FontSize', style.fontSize, ...
        'FontName', style.fontName, ...
        'color', style.backgroundColor, ...
        'XColor', style.fontColor, ...
        'YColor', style.fontColor, ...
        'ZColor', style.fontColor, ...
        'GridColor', style.fontColor);

    axTitle = get(ax, 'title');
    set(axTitle, 'color', style.fontColor);
    
    if ~isempty(hLegend)
        hLegend.Color = style.backgroundColor;
        hLegend.TextColor = style.fontColor;
    end
    
    if ~isempty(hColorbar)
        for i = 1:length(hColorbar)
            hColorbar(i).Color = style.fontColor;
        end
    end

    lines = findobj(gcf,'Type','Line');
    for i = 1:numel(lines)
      set(lines(i), 'LineWidth', style.lineWidth);
    end

    dataObjs = findobj(fig,'-property','YData');
    
    ymin = 0;
    ymax = 1;
    xmin = 0;
    xmax = 1;
    
    for i = 1:length(dataObjs)
        y = get(dataObjs(i), 'YData');
        x = get(dataObjs(i), 'XData');
        
        if i == 1
            ymax = max(max(y));
            ymin = min(min(y));
            xmax = max(max(x));
            xmin = min(min(x));
            
        else
            ymax = max( max(max(y)), ymax);
            ymin = min( min(min(y)), ymin);
            xmax = max( max(max(x)), xmax);
            xmin = min( min(min(x)), xmin);
        end
    end
        ylim([ymin, ymax]);
        xlim([xmin, xmax]);
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
