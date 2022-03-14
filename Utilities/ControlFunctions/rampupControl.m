function val = rampupControl(t, tup, inputI)

    rampupcase = 'sineup';        
    switch rampupcase
        
      case 'sineup'
        val = (t <= tup) .* sineup(0, inputI, 0, tup, t) + (t > tup) .* inputI;
        
      case 'linear'
        val = (t <= tup).*t./tup .* inputI +  (t > tup) .* inputI;
        
    end
    
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
