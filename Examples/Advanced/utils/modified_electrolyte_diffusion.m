function D = modified_electrolyte_diffusion(c, T)
%% modified version of electrolyte conductivity, used to enhance thermal coupling in an example
%% original version is computeDiffusionCoefficient_Ai2020.m
    
    % Calculate diffusion coefficients constant for the diffusion coefficient calculation
    cnst = [ -4.43, -54; 
             -0.22, 0.0 ];

    Tgi = [ 229; 5.0 ];
    
    modified_coef1 = 1; % original is 1

    % Diffusion coefficient, [m^2 s^-1]
    D = 10 .^ ( ( cnst(1,1) + modified_coef1*cnst(1,2) ./ ( T - Tgi(1) - Tgi(2) .* c .* 1e-3) + cnst(2,1) .* ...
                          c .* 1e-3) );

    % D = max(D, 1e-15);
    
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
