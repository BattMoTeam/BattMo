function [] = plotHeatSource(model,states,varargin)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

if model.G.griddim == 2
    xmin = min(model.ThermalModel.G.nodes.coords(:,1));
    xmax = max(model.ThermalModel.G.nodes.coords(:,1));
    ymin = min(model.ThermalModel.G.nodes.coords(:,2));
    ymax = max(model.ThermalModel.G.nodes.coords(:,2));
    ylim([ymin, ymax]);
    xlim([xmin, xmax]);
    
    xlabel('X Position  /  m')
    ylabel('Y Position  /  m')
end

figureHeatSource = plotCellData(model.ThermalModel.G, states.ThermalModel.jHeatSource);
figureHeatSource.EdgeColor = 'none';

end

