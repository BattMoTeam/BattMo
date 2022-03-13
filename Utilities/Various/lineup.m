function [ sw ] = lineup(y1, y2, x1, x2, x)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

dy = y2-y1;
dx = (x2-x1);

m = dy/dx;

x = x-x1;

sw = (x >= x1 & x<= x2) .* m.*x + y1 + ...
    (x>x2) .* y2 + ...
    (x<x1) .* y1;

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
