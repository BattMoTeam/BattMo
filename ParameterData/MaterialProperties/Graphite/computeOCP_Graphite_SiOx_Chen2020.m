function OCP = computeOCP_Graphite_SiOx_Chen2020(sto)
    % LG M50 graphite open circuit potential as a function of stochiometry, fit taken
    % from [1].

    % References
    % ----------
    % .. [1] Chang-Hui Chen, Ferran Brosa Planella, Kieran O’Regan, Dominika Gastol, W.
    % Dhammika Widanage, and Emma Kendrick. "Development of Experimental Techniques for
    % Parameterization of Multi-scale Lithium-ion Battery Models." Journal of the
    % Electrochemical Society 167 (2020): 080534.

    OCP = 1.9793*exp(-39.3631*sto) + 0.2482 - 0.0909*tanh(29.8538*(sto - 0.1234)) - 0.04478*tanh(14.9159*(sto - 0.2769)) - 0.0205*tanh(30.4444*(sto - 0.6103));
    
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
