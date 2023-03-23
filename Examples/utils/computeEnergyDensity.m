function [Emid, Imid, energyDensity, specificEnergy, energy] = computeEnergyDensity(E, I, t, vol, mass, varargin)
% Given time, voltage (E) and current (I), compute the energy that has been produced as a function of time.
% In addition, returns the energyDensity and specificEnergy that are computed using the volume and mass of the cell.
    
    opt = struct('SIoutput', true);
    opt = merge_options(opt, varargin{:});
    
    dt = diff(t);
    Emid = (E(2:end) + E(1:end - 1))/2.0; 
    Imid = (I(2:end) + I(1:end - 1))/2.0; 
    energy = cumsum((Emid.*Imid).*dt); % J

    if opt.SIoutput
        energyDensity = energy/vol; % J/m^3
    else
        energy = energy/hour; % Wh
        energyDensity = energy/vol/1000; %  Wh/L
    end
    
    specificEnergy = energy/mass; % J/kg or Wh/kg
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
