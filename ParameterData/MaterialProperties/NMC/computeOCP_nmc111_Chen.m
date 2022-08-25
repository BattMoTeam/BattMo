function [OCP, dUdT] = computeOCP_nmc111_Chen(c, T, cmax)

    % LG M50 NMC open circuit potential as a function of stochiometry, fit taken
    % from [1].

    % [1] Chang-Hui Chen, Ferran Brosa Planella, Kieran O’Regan, Dominika Gastol, W.
    % Dhammika Widanage, and Emma Kendrick. "Development of Experimental Techniques for
    % Parameterization of Multi-scale Lithium-ion Battery Models." Journal of the
    % Electrochemical Society 167 (2020): 080534.

    sto = c./cmax;
    OCP = -0.8090*sto + 4.4875 - 0.0428*tanh(18.5138*(sto - 0.5542)) - 17.7326*tanh(15.7890*(sto - 0.3117)) + 17.5842*tanh(15.9308*(sto - 0.3120));
    dUdT = 0;
    
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
