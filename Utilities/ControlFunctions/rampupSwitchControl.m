function [val, ctrltype] = rampupSwitchControl(t, tup, I, E, inputI, inputE)
% CURRENTSOURCE Summary of this function goes here
%   Detailed explanation goes here

    if(inputI<0)
        %notuseVolt = ((E - inputE)<1e-3);
        %notuseVolt = not(useVolt);
        useCurrent = (E < inputE);
    else
        %notuseVolt = ((E - inputE)>-1e-3);
        useCurrent = (E > inputE);
        %useVolt = not(useVolt);
    end
    %if not(useVolt) || (t < tup)
    if  useCurrent | (t < tup)   
        % We control with current
        ctrltype = 'I';
        
        rampupcase = 'sineup';        
        switch rampupcase
            
          case 'sineup'
            val = (t <= tup) .* sineup(0, inputI, 0, tup, t) + (t > tup) .* inputI;
            
          case 'linear'
            val = (t <= tup).*t./tup .* inputI +  (t > tup ) .* inputI;
            
        end
    
    else
        
        % We control with current
        ctrltype = 'E';

        val = inputE;
        
    end
    
end



%{
Copyright 2009-2021 SINTEF Industry, Sustainable Energy Technology
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
