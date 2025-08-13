function dUdT = computeEntropyChange_NMC111(theta)
    
    % Calculate the entropy change at the given lithiation
    
    coeff1 = [0.199521039        , ...
                   - 0.928373822      , ...
                   + 1.364550689000003, ...
                   - 0.611544893999998];
    
    coeff2 = [1                  , ...
                   - 5.661479886999997, ...
                   + 11.47636191      , ... 
                   - 9.82431213599998 , ...
                   + 3.048755063];
    
    dUdT = -1e-3.*polyval(coeff1(end:-1:1),theta)./ polyval(coeff2(end:-1:1), theta);
    
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
