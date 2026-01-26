function D = computeDiffusionCoefficient_Ai2020(c, T)
% Diffusivity of LiPF6 in EC:DMC as a function of ion concentration.
% 
% References
% ----------
% .. [1] Ai, W., Kraft, L., Sturm, J., Jossen, A., & Wu, B. (2020).
% Electrochemical Thermal-Mechanical Modelling of Stress Inhomogeneity
% in Lithium-Ion Pouch Cells. Journal of The Electrochemical Society,
% 167(1), 013512. DOI: 10.1149/2.0122001JES.
%
% Calculate diffusion coefficients constant for the diffusion coefficient calculation
    
    cnst = [ -4.43, -54; 
             -0.22, 0.0 ];

    Tgi = [ 229; 5.0 ];
    
    % Diffusion coefficient, [m^2 s^-1]
    D = 1e-4 .* 10 .^ ( ( cnst(1,1) + cnst(1,2) ./ ( T - Tgi(1) - Tgi(2) .* c .* 1e-3) + cnst(2,1) .* ...
                          c .* 1e-3) );
    
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
