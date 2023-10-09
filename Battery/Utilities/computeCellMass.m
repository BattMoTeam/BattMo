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

        switch model.(elde).(co).activematerial_type
            
          case 'default'
            
            rho  = model.(elde).(co).density;
            vols = model.(elde).(co).G.cells.volumes;
            frac = model.(elde).(co).volumeFraction;
            
            masses.(elde).(co).val  = sum(rho.*vols.*frac);
            volumes.(elde).(co).val = sum(vols.*frac);
            
          case 'composite'
            
            gr = 'FirstMaterial';
            si = 'SecondMaterial';

            mats = {gr, si};
            
            masses.(elde).(am).val  = 0;
            volumes.(elde).(am).val = 0;
            
            for imat = 1 : numel(mats)

                mat = mats{imat};
                rho  = model.(elde).(am).(mat).(itf).density;
                vols = model.(elde).(am).(mat).G.cells.volumes;
                vf = model.(elde).(am).volumeFraction;
                amvf = model.(elde).(am).(mat).activeMaterialFraction;
                
                masses.(elde).(am).(mat).val = sum(rho.*vols.*vf.*amvf);
                masses.(elde).(am).val = masses.(elde).(am).val + masses.(elde).(am).(mat).val;

                volumes.(elde).(am).(mat).val = sum(vols.*vf.*amvf);
                volumes.(elde).(am).val = volumes.(elde).(am).val + volumes.(elde).(am).(mat).val;

            end
            
          otherwise
            
            error('electrode_case not recognized');
            
        end
        
        mass = mass + masses.(elde).(co).val;
        
        if model.include_current_collectors

            rho  = model.(elde).(cc).density;
            vols = model.(elde).(cc).G.cells.volumes;
            
            masses.(elde).(cc).val  = sum(rho.*vols);
            volumes.(elde).(cc).val = sum(vols);

            mass = mass + masses.(elde).(cc).val;
            
        end
        
    end
    
    rho  = model.(elyte).density;
    vols = model.(elyte).G.cells.volumes;
    frac = model.(elyte).volumeFraction;
    
    masses.(elyte).val  = sum(rho.*vols.*frac);
    volumes.(elyte).val = sum(vols.*frac);
    
    mass = mass + masses.(elyte).val;
    
    rho  = model.(sep).density;
    vols = model.(sep).G.cells.volumes;
    frac = (1 - model.(sep).porosity);
    
    masses.(sep).val  = sum(rho.*vols.*frac);
    volumes.(sep).val = sum(vols.*frac);

    mass = mass + masses.(sep).val;

    mass = mass + opt.packingMass;

    volumes.val = sum(model.G.cells.volumes) + opt.packingVolume;
    
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
