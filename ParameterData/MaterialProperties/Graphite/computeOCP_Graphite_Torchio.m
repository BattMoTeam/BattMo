function OCP = computeOCP_Graphite_Torchio(theta)
%   Calculate the equilibrium open cirucuit potential of
%   graphite according to the model used by Torchio et al [1].

    % Calculate the open-circuit potential at the reference temperature for the given lithiation
    OCP = (0.7222 ...
           + 0.1387 .* theta ...
           + 0.0290 .* theta.^0.5 ...
           - 0.0172 ./ theta ... 
           + 0.0019 ./ theta.^1.5 ...
           + 0.2808 .* exp(0.9 - 15.*theta) ... 
           - 0.7984 .* exp(0.4465.*theta - 0.4108));

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
