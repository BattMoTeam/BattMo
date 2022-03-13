function windingnumbers = computeWindingNumbers(fractions, rInner, w, nwindings)

    x  = linspace(0, 2*pi*nwindings*rInner, 10000);
    dx = diff(x);

    params.w = w;
    params.rInner = rInner;
    dlperdx = spiralDeltaLength(x, 0, params);
    l = cumsum(dlperdx(1 : (end - 1)).*dx);
    
    totallength = l(end);
    
    for ind = 1 : numel(fractions)
        xind = find(l/totallength > fractions(ind), 1, 'first');
        windingnumbers(ind) = ceil(x(xind)/(2*pi*rInner));
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
along with BattMog.  If not, see <http://www.gnu.org/licenses/>.
%}
