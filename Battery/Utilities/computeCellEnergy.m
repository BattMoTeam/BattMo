function [energy, extras] = computeCellEnergy(model, varargin)
%% Compute the cell energy.
%
% If no DRate is given, it corresponds to the maximum theoritical cell energy, at infinitly small DRate. Otherwise,
% simulate the discharge and gives the corresponding energy.
%
% The extras structure provides more detailed information with the fields
% - energy
% - dischargeFunction % function handler giving the voltage as a function of the state of charge
% - E                 % Voltage output (raw computation data in case DRate is given)
% - I                 % Current output (raw computation data in case DRate is given)
% - time              % time output (raw computation data in case DRate is given)
    
    opt = struct('capacities' , [] , ...
                 'temperature', 298, ...
                 'DRate', []);
    [opt, extra] = merge_options(opt, varargin{:});

    ne  = 'NegativeElectrode';
    pe  = 'PositiveElectrode';
    am  = 'ActiveMaterial';
    co  = 'Coating';
    itf = 'Interface';
    sd  = 'SolidDiffusion';

    if isempty(opt.DRate)

        if model.(ne).(co).activeMaterialModelSetup.composite || model.(pe).(co).activeMaterialModelSetup.composite

            [energy, extras] = advancedComputeCellEnergy(model);
            return
            
        end
        
        if isempty(opt.capacities)
            [~, capacities] = computeCellCapacity(model);
        else
            capacities = opt.capacities
        end

        T = opt.temperature;

        capacity = min(capacities.(ne), capacities.(pe));


        th100 = 'guestStoichiometry100';
        th0   = 'guestStoichiometry0';
        sc    = 'saturationConcentration';

        N = 1000;

        eldes = {ne, pe};
        
        for ielde = 1 : numel(eldes)

            elde = eldes{ielde};
            
            smax = capacity./capacities.(elde);
            
            itfmodel = model.(elde).(co).(am).(itf);
            c0 = itfmodel.(th100)*itfmodel.(sc);
            cT = itfmodel.(th0)*itfmodel.(sc);            
            
            s = smax.*linspace(0, 1, N + 1)';

            c    = (1 - s).*c0 + s.*cT;
            cmax = model.(elde).(co).(am).(itf).(sc);
            
            f = model.(elde).(co).(am).(itf).computeOCP(c(1 : end - 1)/cmax);
            
            % function handler
            fs{ielde} = @(s) model.(elde).(co).(am).(itf).computeOCP(((1 - s).*c0 + s.*cT)/cmax);
            
            energies{ielde} = capacities.(elde)*smax/N*sum(f);
            
        end
        
        energy = (energies{2} - energies{1});

        dischargeFunction = @(s) (fs{2}(s) - fs{1}(s));

        extras = struct('energy'           , energy, ...
                        'dischargeFunction', dischargeFunction);
        
    else

        DRate  = opt.DRate;
        extras = computeCellEnergyGivenDrate(model, DRate, extra{:});
        
    end

    energy = extras.energy;
    
end



%{
Copyright 2021-2024 SINTEF Industry, Sustainable Energy Technology
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
