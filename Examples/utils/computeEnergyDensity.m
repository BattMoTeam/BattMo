function [Emid, Tmid, energyDensity, specificEnergy, Energy] = computeEnergyDensity(E, I, T, t, vol, mass)
    dt = diff(t);
    Emid = (E(2:end) + E(1:end - 1))/2.0; 
    Tmid = (T(2:end) + T(1:end - 1))/2.0;
    Imid = (I(2:end) + I(1:end - 1))/2.0; 
    Energy = cumsum((Emid.*Imid).*dt);
    Energy = Energy/hour; % W h
    energyDensity = Energy/vol/1000; %  W h/L
    specificEnergy = Energy/mass; %  W h/kg
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
