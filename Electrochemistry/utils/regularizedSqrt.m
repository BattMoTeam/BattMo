function y = regularizedSqrt(x, th)

%
%
% SYNOPSIS:
%   function y = regularizedSqrt(x, th)
%
% DESCRIPTION: returns regularized square root by using linear interporation between (0, 0) and (th, sqrt(th))
%
% PARAMETERS:
%   x  - input values
%   th - threshold
%
% RETURNS:
%   y - output values
%

    y = x; % quick way to create y of same dimension as x and also preserved AD
    ind = (x <= th);
    y(~ind) = (x(~ind)).^0.5;
    y(ind) = x(ind)/th*sqrt(th);
    
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
