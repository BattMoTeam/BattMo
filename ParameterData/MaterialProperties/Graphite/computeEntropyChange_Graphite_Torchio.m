function dUdT = computeEntropyChange_Graphite_Torchio(theta)
% Graphite entropy change from the fit in Table II of Torchio et al. [1].
%
% References
% ----------
% .. [1] Torchio, M., Magni, L., Gopaluni, R. B., Braatz, R. D., and Raimondo, D. M. (2016).
%    LIONSIMBA: A Matlab Framework Based on a Finite Volume Model Suitable for Li-Ion Battery
%    Design, Simulation, and Control. Journal of The Electrochemical Society, 163(7),
%    A1192-A1205. DOI: 10.1149/2.0291607jes.
    
    coeff1 = [0.005269056 ,...
              + 3.299265709,...
              - 91.79325798,...
              + 1004.911008,...
              - 5812.278127,...
              + 19329.75490,...
              - 37147.89470,...
              + 38379.18127,...
              - 16515.05308];
    
    coeff2= [1, ...
             - 48.09287227,...
             + 1017.234804,...
             - 10481.80419,...
             + 59431.30000,...
             - 195881.6488,...
             + 374577.3152,...
             - 385821.1607,...
             + 165705.8597];
    
    dUdT = 1e-3.*polyval(coeff1(end:-1:1),theta)./ polyval(coeff2(end:-1:1),theta);
    
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
