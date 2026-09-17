function D = computeDiffusionCoefficient_Chen2020(c, T)
% Electrolyte transport fit used with the Chen et al. parameter set [1]; the underlying
% electrolyte measurements are from Nyman et al. [2].
%
% References
% ----------
% .. [1] Chen, C.-H., Brosa Planella, F., O’Regan, K., Gastol, D., Widanage, W. D., and Kendrick,
%    E. (2020). Development of Experimental Techniques for Parameterization of Multi-scale
%    Lithium-ion Battery Models. Journal of The Electrochemical Society, 167, 080534. DOI:
%    10.1149/1945-7111/ab9050.
% .. [2] Nyman, A., Behm, M., and Lindbergh, G. (2008). Electrochemical characterisation and
%    modelling of the mass transport phenomena in LiPF6-EC-EMC electrolyte. Electrochimica
%    Acta, 53(22), 6356-6365. DOI: 10.1016/j.electacta.2008.04.023.
    
    c = c./1000;
    D = 8.794 .*10^(-11) .*c .^2 - 3.972 .*10^(-10) .*c + 4.862*10^(-10);

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
