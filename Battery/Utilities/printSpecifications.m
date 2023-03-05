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
