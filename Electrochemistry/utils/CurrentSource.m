function res = CurrentSource(t, tup, tlim, jnom)
% CURRENTSOURCE Summary of this function goes here
%   Detailed explanation goes here

N = length(jnom);

if N == 1
    
    rampupcase = 'sineup';

    switch rampupcase
      
      case 'sineup'
        res = (t <= tup) .* sineup(0, jnom, 0, tup, t) + (t > tup && t <= tlim) .* jnom;
      
      case 'linear'
        res = (t <= tup).*t./tup .* jnom +  (t > tup && t <= tlim) .* jnom;
        
    end
    
else
    
    tboundslow = zeros(N, 1);
    tboundsup = zeros(N, 1);
    tboundsup(1) = tlim(1);
    
    for i = 2:N
        tboundslow(i) = tboundslow(i-1) + tlim(i-1);
        tboundsup(i) = tboundsup(i-1) + tlim(i);
    end
    
    if t == tboundsup(end)
        step = N;
    else
    step = find(t >= tboundslow & t < tboundsup);
    end
    
    shift = tboundslow(step);
    if step == 1
        res = (t <= tup(step)) .* sineup(0, jnom(step), 0, tup(step), t) + ...
            ... Initial Ramp Up
            (t > tup(step) & t <= tlim(step)) .* jnom(step);
    else
        res = (t>= shift & t <= (tup(step) + shift) ) .* sineup(jnom(step-1), jnom(step), tboundslow(step), tboundslow(step) + tup(step), t) + ...
            ... Initial Ramp Up
            (t > tboundslow(step) + tup(step) & t <= tboundsup(step)) .* jnom(step);
    end
end
        

% + ...
%     ... Hold nominal current
%     (t>ZAB.td0 && t<= (ZAB.td0 + ZAB.dts)) .* sineup (ZAB.jd, -ZAB.jc, ZAB.td0, (ZAB.td0 + ZAB.dts), t) + ...
%     ... Cycle 1 Ramp to Charge
%     (t>(ZAB.td0 + ZAB.dts) && t <= (ZAB.td0 + ZAB.tc)) .* (-1) .* ZAB.jc + ...
%     ... Cycle 1 Hold Charge
%     (t>(ZAB.td0 + ZAB.tc) && t <= (ZAB.td0 + ZAB.tc + ZAB.dts)) .* sineup(-ZAB.jc, ZAB.jd, (ZAB.td0 + ZAB.tc), (ZAB.td0 + ZAB.tc + ZAB.dts), t) + ...
%     ... Cycle 1 Ramp to Discharge
%     (t>(ZAB.td0 + ZAB.tc + ZAB.dts) && t<= (ZAB.td0 + ZAB.tc + ZAB.td)) .* ZAB.jd;
%     ... Cycle 1 Hold Discharge

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
