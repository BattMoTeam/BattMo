function [Z_real, Z_imag] = load_nyquist(params, omega)

    R_0 = params(1);
    R_1 = params(2);
    C_1 = params(3);
    R_2 = params(4);
    C_2 = params(5);
    
    Z_real = R_0 + R_1 ./ (1+(R_1* C_1.*omega).^2) +   R_2 ./ (1+(R_2* C_2.*omega).^2);
    Z_imag = -(R_1*R_1*C_1.*omega ./ (1+(R_1* C_1.*omega).^2))  ...
               - (R_2*R_2*C_2.*omega ./ (1+(R_2* C_2.*omega).^2));

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
