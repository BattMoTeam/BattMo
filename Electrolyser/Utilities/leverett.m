function J = leverett(par, s)
    %LEVERETT calculates the leverett J-function of a porous
    %domain. 
    %   obj.leverett() calculates the leverett J-function
    %   for a liquid saturation of s, considering the properties of
    %   the porous medium.
    
    J0  = par(1);
    A1  = par(2);
    A2  = par(3);
    B1  = par(4);
    B2  = par(5);

    J   =   J0 + ...
            A1.*exp(B1.*(s - 0.5)) - ...
            A2.*exp(B2.*(s - 0.5));
                      
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
