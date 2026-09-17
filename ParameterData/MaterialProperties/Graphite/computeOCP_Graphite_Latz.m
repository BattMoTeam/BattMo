function OCP = computeOCP_Graphite_Latz(theta)
% computeOCP_Graphite_Latz
%   Computes the equilibrium open-circuit potential (OCP) of
%   graphite based on the model from Hein, Danner, and Latz [1].
%
%   Inputs:
%       theta - state of charge
%    
%   Output:
%       OCP   - open-circuit potential at 298.15 K [V vs Li/Li+]
%
% References
% ----------
% .. [1] Hein, S., Danner, T., and Latz, A. (2020). An Electrochemical Model of Lithium Plating and
%    Stripping in Lithium Ion Batteries. ACS Applied Energy Materials, 3(9), 8519-8531. DOI:
%    10.1021/acsaem.0c01155.

    % Modified hyperbolic tangent function used in the original paper
    tanhmod = @(x) (exp(20 .* x) - exp(-x)) ./ (exp(-x) + exp(x));

    % OCP expression at reference temperature (T = 298.15 K)
    OCP = 0.6379 ...
        + 0.5416 .* exp(-305.5309 .* theta) ...
        + 0.044  .* tanh((-theta - 0.1958) ./ 0.1088) ...
        - 0.1978 .* tanhmod((theta - 0.99) ./ 0.05) ...
        - 0.6875 .* tanh((theta + 0.0117) ./ 0.0529) ...
        - 0.0175 .* tanh((theta - 0.5692) ./ 0.0875);


end

%{
Copyright 2021-2026 SINTEF Industry, Sustainable Energy Technology
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
