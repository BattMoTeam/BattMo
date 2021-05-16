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

