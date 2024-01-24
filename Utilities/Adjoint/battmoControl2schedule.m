function schedule = battmoControl2schedule(u, schedule, scaling)
    %% Convert control vector u to schedule 
    
    nc = numel(schedule.control);
    [umin, umax] = deal(scaling.boxLims(:,1), scaling.boxLims(:,2));
    c = 0;
    
    for cs = 1 : nc

        c       = c + 1;
        cutoffV = 3.0;
        tup     = 0.1;

        Imax =  u(c) *(umax - umin) + umin;

        schedule.control(cs).src = @(time, I, E) rampupSwitchControl(time, tup, I, E, Imax, cutoffV);
        
    end
    
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
