function [OCP, dUdT] = computeOCP_LFP_Gerver2011(c, T, cmax)
% LFP open-circuit potential. The filename attributes the fit to Gerver (2011); [1] is a related
% publication. The exact source of these coefficients has not been verified.
%
% References
% ----------
% .. [1] Gerver, R. E., and Meyers, J. P. (2011). Three-Dimensional Modeling of Electrochemical
%    Performance and Heat Generation of Lithium-Ion Batteries in Tabbed Planar Configurations.
%    Journal of The Electrochemical Society, 158(7), A835-A843. DOI: 10.1149/1.3591799.

    error("function is not compatible with new function interface");
    
    OCP = 3.41285712e+00 ...
          - 1.49721852e-02 * c/cmax ...
          + 3.54866018e+14 * exp(-3.95729493e+02 * c/cmax) ...
          - 1.45998465e+00 * exp(-1.10108622e+02 * (1 - c/cmax));
    
    dUdT = 0;
    
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
