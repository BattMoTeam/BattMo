function [] = plotTemperature(model, states, varargin)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

if model.G.griddim == 2
    xmin = min(model.ThermalModel.G.nodes.coords(:,1));
    xmax = max(model.ThermalModel.G.nodes.coords(:,1));
    ymin = min(model.ThermalModel.G.nodes.coords(:,2));
    ymax = max(model.ThermalModel.G.nodes.coords(:,2));
    ylim([ymin, ymax]);
    xlim([xmin, xmax]);
    
    xlabel(gca, 'X Position  /  m')
    ylabel(gca, 'Y Position  /  m')
end

figureTemperature = plotCellData(model.ThermalModel.G, states.ThermalModel.T);
figureTemperature.EdgeColor = 'none';
end




%{
Copyright 2009-2021 SINTEF Industry, Sustainable Energy Technology
and SINTEF Digital, Mathematics & Cybernetics.

This file is part of The Battery Modeling Toolbox BatMo

BatMo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BatMo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BatMo.  If not, see <http://www.gnu.org/licenses/>.
%}
