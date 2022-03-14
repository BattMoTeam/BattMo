classdef FigureStyle
    %PLOTSTYLE Defines the style properties of a figure window
    %   Detailed explanation goes here
    
    properties      
        width                       % Canvas width
        height                      % Canvas height
        fontName                    % Font name
        fontSize                    % Font size
        fontColor                   % Font color
        backgroundColor             % Figure background color
        colormap_temperature        % Colormap for temperature values
        colormap_heatSource         % Colormap for heat source values
        colormap_concentration      % Colormap for concentration values
        colormap_reactionRate       % Colormap for reactionrate values
        colormap_electricPotential  % Colormap for electric potential values
        colorOrder                  % Color order for discrete datasets
    end
    
    methods
        function obj = FigureStyle(varargin)
            %FigureStyle Construct an instance of this class
            %   Detailed explanation goes here
            
            %% Parse inputs            
            defaultTheme        = 'light';
            expectedTheme       = {'dark', 'light', 'blue'};
            defaultSize         = 'A4';
            expectedSize        = {'A4', '10x10'};
            defaultOrientation  = 'landscape';
            expectedOrientation = {'landscape', 'portrait'};

            p = inputParser;
            addParameter(p, 'theme', defaultTheme, ...
                @(x) any(validatestring(x,expectedTheme)));
            addParameter(p, 'size', defaultSize, ...
                @(x) any(validatestring(x,expectedSize)));
            addParameter(p, 'orientation', defaultOrientation, ...
                @(x) any(validatestring(x,expectedOrientation)));
            parse(p, varargin{:});
            
            %% Define test properties
            obj.fontName = 'Helvetica';
            obj.fontSize = 24;
            
            %% Define size properties
            if strcmpi(p.Results.size, 'A4')
                x1 = 21;    % cm
                x2 = 29.7;  % cm
                if strcmpi(p.Results.orientation, 'landscape')
                    obj.width   = max(x1, x2);
                    obj.height  = min(x1, x2);
                elseif strcmpi(p.Results.orientation, 'portrait')
                    obj.width   = min(x1, x2);
                    obj.height  = max(x1, x2);
                end
            end
            
            %% Define color properties for standard plot themes
            % Default values
            obj.colormap_temperature        = 'parula';
            obj.colormap_heatSource         = 'parula';
            obj.colormap_concentration      = 'parula';
            obj.colormap_reactionRate       = 'parula';
            obj.colormap_electricPotential  = 'parula';
            
            if strcmpi(p.Results.theme,'dark')
                obj.backgroundColor         = coal();
                obj.fontColor               = 'w';
                obj.colormap_temperature    = 'inferno';
            elseif strcmpi(p.Results.theme,'light')
                obj.backgroundColor         = 'w';
                obj.fontColor               = 'k';
                obj.colormap_temperature    = 'inferno';
            elseif strcmpi(p.Results.theme,'blue')
                obj.backgroundColor         = navyblue();
                obj.fontColor               = 'w';
                obj.colormap_temperature    = 'inferno';
            end
            
        end
        
        function [] = applyFigureStyle(obj, H)
           
            set(H, ...
                'units','centimeter','position',[1, 1, obj.width, obj.height,], ...
                'color', obj.backgroundColor)
            
            set(gca, ...
                'FontSize',obj.fontSize, ...
                'FontName', obj.fontName, ...
                'color', obj.backgroundColor, ...
                'XColor', obj.fontColor, ...
                'YColor', obj.fontColor, ...
                'GridColor', obj.fontColor);
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
