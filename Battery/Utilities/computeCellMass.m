function [mass, masses, volumes] = computeCellMass(model, varargin)

    opt = struct('packingMass', 0, ...
                 'packingVolume', 0);
    opt = merge_options(opt, varargin{:});
    
    
    elyte = 'Electrolyte';
    sep   = 'Separator';
    ne    = 'NegativeElectrode';
    pe    = 'PositiveElectrode';
    co    = 'Coating';
    am    = 'ActiveMaterial';
    itf   = 'Interface';
    cc    = 'CurrentCollector';
    
    eldes = {ne, pe};
    
    mass = 0; % total mass
    
    % We first compute the mass of the electrodes
    for ind = 1 : numel(eldes)
        
        elde = eldes{ind};

        rho  = model.(elde).(co).effectiveDensity;
        vols = model.(elde).(co).G.getVolumes();
        frac = model.(elde).(co).volumeFraction;
        
        masses.(elde).(co).val  = sum(rho.*vols);
        volumes.(elde).(co).val = sum(vols.*frac);
        
        mass = mass + masses.(elde).(co).val;
        
        if model.include_current_collectors

            rho  = model.(elde).(cc).density;
            vols = model.(elde).(cc).G.getVolumes();
            
            masses.(elde).(cc).val  = sum(rho.*vols);
            volumes.(elde).(cc).val = sum(vols);

            mass = mass + masses.(elde).(cc).val;
            
        end
        
    end
    
    rho  = model.(elyte).density;
    vols = model.(elyte).G.getVolumes();
    frac = model.(elyte).volumeFraction;
    
    masses.(elyte).val  = sum(rho.*vols.*frac);
    volumes.(elyte).val = sum(vols.*frac);
    
    mass = mass + masses.(elyte).val;
    
    rho  = model.(sep).density;
    vols = model.(sep).G.getVolumes();
    frac = (1 - model.(sep).porosity);
    
    masses.(sep).val  = sum(rho.*vols.*frac);
    volumes.(sep).val = sum(vols.*frac);

    mass = mass + masses.(sep).val;

    mass = mass + opt.packingMass;

    volumes.val = sum(model.G.getVolumes()) + opt.packingVolume;
    
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
