function D = computeDiffusionCoefficient_Nyman2008(c, T)
% References
% ----------
% .. [1] Nyman, A., Behm, M., and Lindbergh, G. (2008). Electrochemical characterisation and
%    modelling of the mass transport phenomena in LiPF6-EC-EMC electrolyte. Electrochimica
%    Acta, 53(22), 6356-6365. DOI: 10.1016/j.electacta.2008.04.023.

    D = 8.794e-11*(c/1000).^2 - 3.972e-10*(c/1000) + 4.862e-10;

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
