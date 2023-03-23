function printSpecifications(model, varargin)

    opt = struct('packingMass', 0);
    opt = merge_options(opt, varargin{:});
    
    [mass, masses] = computeCellMass(model, 'packingMass', opt.packingMass); % mass in kg

    G = model.G;
    vol = sum(G.cells.volumes); % volume in cubic meter

    [cap, cap_neg, cap_pos, specificEnergy] = computeCellCapacity(model, 'packingMass', opt.packingMass); 
    % cap, cap_neg, cap_pos in Coulomb
    % specificEnergy in Joule/kg

    specificEnergy = specificEnergy/hour;

    % electrode capacities
    ne      = 'NegativeElectrode';
    pe      = 'PositiveElectrode';
    am      = 'ActiveMaterial';
    itf     = 'Interface';
    
    fprintf('%16s: %g [g]\n','packing mass', opt.packingMass/gram);
    fprintf('%16s: %g [kg]\n', 'mass', mass);
    fprintf('%16s: %g [l]\n', 'volume', vol/(1*litre));
    fprintf('%16s: %g [Ah]\n', 'Capacity', cap/(1*hour));
    % fprintf('%16s: %g [Ah]\n', 'Capacity (neg)', cap_neg/(1*hour));
    % fprintf('%16s: %g [Ah]\n', 'Capacity (pos)', cap_pos/(1*hour));
    fprintf('%16s: %g [Wh]\n', 'Energy', specificEnergy*mass);
    fprintf('%16s: %g [Wh/l]\n', 'Energy Density', ((specificEnergy*mass)/(vol/(1*litre))));
    fprintf('%16s: %g [Wh/kg]\n', 'Specific Energy', specificEnergy);
    
end



%{
Copyright 2021-2023 SINTEF Industry, Sustainable Energy Technology
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
