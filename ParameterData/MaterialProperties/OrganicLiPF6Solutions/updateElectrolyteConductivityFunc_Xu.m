function conductivity = updateElectrolyteConductivityFunc_Xu(c, T)
    
    conductivityFactor = 1e-4;
    
    %cnst = [-10.5   , 0.074    , -6.96e-5; ...
    %        0.668e-3, -1.78e-5 , 2.80e-8; ...
    %        0.494e-6, -8.86e-10, 0];            
            
    
    % Ionic conductivity, [S m^-1]
    %conductivity = conductivityFactor.* c .*( polyval(cnst(end:-1:1,1),c) + polyval(cnst(end:-1:1,2),c) .* T + ...
    %                                          polyval(cnst(end:-1:1,3),c) .* T.^2).^2;
    %From cideMOD
    conductivity=c.*1e-4*1.2544.* ...
        (-8.2488+0.053248.*T-2.987e-5*(T.^2)+ ...
        0.26235e-3.*c-9.3063e-6.*c.*T+ ...
        8.069e-9.*c.*T.^2+ ...
        2.2002e-7.*c.^2-1.765e-10.*T.*c.^2);
    
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
