function [] = plotElectricPotential(model, states, varargin)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

%% Parse inputs
defaultDomain       =   'Electrolyte';
expectedDomain      = { 'NegativeActiveMaterial', ...
                        'NegativeElectrodeCurrentCollector', ...
                        'PositiveActiveMaterial', ...
                        'PositiveElectrodeCurrentCollector', ...
                        'Electrolyte' };

p = inputParser;
addRequired(p, 'model')
addRequired(p, 'states')
addParameter(p, 'domain', defaultDomain, ...
    @(x) any(validatestring(x,expectedDomain)));
parse(p, model, states, varargin{:});

if strcmpi(p.Results.domain, 'Electrolyte')
    domain = model.Electrolyte;
    state = states.Electrolyte;
elseif strcmpi(p.Results.domain, 'NegativeActiveMaterial')
    domain = model.NegativeElectrode.ActiveMaterial;
    state = states.NegativeElectrode.ActiveMaterial;
elseif strcmpi(p.Results.domain, 'NegativeElectrodeCurrentCollector')
    domain = model.NegativeElectrode.CurrentCollector;
    state = states.NegativeElectrode.CurrentCollector;
elseif strcmpi(p.Results.domain, 'PositiveActiveMaterial')
    domain = model.PositiveElectrode.ActiveMaterial;
    state = states.PositiveElectrode.ActiveMaterial;
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
