function [] = plotElectricPotential(model, states, varargin)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

%% Parse inputs
defaultDomain       =   'Electrolyte';
expectedDomain      = { 'NegativeElectrodeActiveComponent', ...
                        'NegativeElectrodeCurrentCollector', ...
                        'PositiveElectrodeActiveComponent', ...
                        'PositiveElectrodeCurrentCollector', ...
                        'Electrolyte' };

p = inputParser;
addRequired(p, 'model')
addRequired(p, 'states')
addOptional(p, 'domain', defaultDomain, ...
    @(x) any(validatestring(x,expectedDomain)));
parse(p, model, states, varargin{:});

if strcmpi(p.Results.domain, 'Electrolyte')
    domain = model.Electrolyte;
    state = states.Electrolyte;
elseif strcmpi(p.Results.domain, 'NegativeElectrodeActiveComponent')
    domain = model.NegativeElectrode.ElectrodeActiveComponent;
    state = states.NegativeElectrode.ElectrodeActiveComponent;
elseif strcmpi(p.Results.domain, 'NegativeElectrodeCurrentCollector')
    domain = model.NegativeElectrode.CurrentCollector;
    state = states.NegativeElectrode.CurrentCollector;
elseif strcmpi(p.Results.domain, 'PositiveElectrodeActiveComponent')
    domain = model.PositiveElectrode.ElectrodeActiveComponent;
    state = states.PositiveElectrode.ElectrodeActiveComponent;
elseif strcmpi(p.Results.domain, 'PositiveElectrodeCurrentCollector')
    domain = model.PositiveElectrode.CurrentCollector;
    state = states.PositiveElectrode.CurrentCollector;
end


if model.G.griddim == 2
    xmin = min(domain.G.nodes.coords(:,1));
    xmax = max(domain.G.nodes.coords(:,1));
    ymin = min(domain.G.nodes.coords(:,2));
    ymax = max(domain.G.nodes.coords(:,2));
    ylim([ymin, ymax]);
    xlim([xmin, xmax]);
    
    xlabel(gca, 'X Position  /  m')
    ylabel(gca, 'Y Position  /  m')
end

plotHandle = plotCellData(domain.G, state.phi);
plotHandle.EdgeColor = 'none';
end
