function dUdT = computeEntropyChange_LCO(theta)

    % This data is taken from PyBamm. 

    % Lithium Cobalt Oxide (LiCO2) entropic change in open-circuit potential (OCP) at
    % a temperature of 298.15K as a function of the stoichiometry. The fit is taken
    % from Scott Moura's FastDFN code [1].

    % References
    % ----------
    % https://github.com/scott-moura/fastDFN

    
    stretch = 1.062;
    theta = stretch * theta;
    % Original parametrization was expressed in terms of c_s_max, but we want to
    % express it in terms of thetaichiometry only
    c_s_max = 51217.9257309275;
    
    dUdT = 0.07645 * (-54.4806 / c_s_max) .* ((1.0 ./ cosh(30.834 - 54.4806 * theta)) .^ 2) ...
           + 2.1581 * (-50.294 / c_s_max) .* ((cosh(52.294 - 50.294 * theta)) .^ (-2)) ...
           + 0.14169 * (19.854 / c_s_max) .* ((cosh(11.0923 - 19.8543 * theta)) .^ (-2)) ...
           - 0.2051 * (5.4888 / c_s_max) .* ((cosh(1.4684 - 5.4888 * theta)) .^ (-2)) ...
           - (0.2531 / 0.1316 / c_s_max) .* ((cosh((-theta + 0.56478) / 0.1316)) .^ (-2)) ...
           - (0.02167 / 0.006 / c_s_max) .* ((cosh((theta - 0.525) / 0.006)) .^ (-2));

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
