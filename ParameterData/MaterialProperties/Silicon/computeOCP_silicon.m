function [OCP, dUdT] = computeOCP_silicon(c, T, cmax)
% Calculate the equilibrium open cirucuit potential of silicon according to the model
% used by Li Haoliang, Bo Lu, Yicheng Song and Junqian Zhang 2017

% The thermal model has not been modified, the values for coeff1 and
% coeff2 are the one for graphite

% It is a swelling Material --> theta cannot just be expressed by c/cmax,
% necessary to come back to the real definition in term of matter
% quantities : N/Nmax; Here, it is only ok for delithiation.

    
    Tref = 298.15;  % [K]

    %reporting useful values from the silicon json file (schould be
    %improved by directly linking the values)

    theta100 = 0.89943252232;
    theta0 = 0.03212259;
    radius0 = 60e-09;
    MolarMassSi   = 28.0855e-3;
    densitySi = 2330;
    molarVolumeLi = 8.8 * 1E-6;
    pi = 3.141592653589793 ;
    
    molarVolumeSi = MolarMassSi/densitySi;
    % same calcul as done in updateRadius_delithiation in the class
    % SwellingMaterial.
    c_ratio = c./(theta100 .*cmax);

    ratio_delith = (3.75.*molarVolumeLi)./(molarVolumeSi+3.75.*molarVolumeLi);
    radius = radius0 .* ((1-ratio_delith) ./ ((1./3.8) + ratio_delith.*c_ratio)) .^ (1/3);

    N = c .* 4/3 .* pi .* radius .^3;
    Nmax = cmax .* pi .* radius0 .^3;

    theta = N./Nmax;


    z = (theta - theta0)/(theta100 - theta0);
    %capping to logical values
    z = min(z, 1);
    z = max(z, 0);
    
    
    % Calculate the open-circuit potential at the reference temperature for the given lithiation
    refOCP = (0.62 ...
              - 1.94 .* z ...
              + 5.8 .*  z.^2 ...
              - 7.13 .* z.^3 ... 
              - 1.8 .*  z.^4 ...
              + 9.34 .* z.^5 ... 
              - 4.76 .* z.^6);
   
    coeff1 = [0.005269056 ,...
              + 3.299265709,...
              - 91.79325798,...
              + 1004.911008,...
              - 5812.278127,...
              + 19329.75490,...
              - 37147.89470,...
              + 38379.18127,...
              - 16515.05308];
    
    coeff2= [1, ...
             - 48.09287227,...
             + 1017.234804,...
             - 10481.80419,...
             + 59431.30000,...
             - 195881.6488,...
             + 374577.3152,...
             - 385821.1607,...
             + 165705.8597];
    
    dUdT = 1e-3.*polyval(coeff1(end:-1:1),theta)./ polyval(coeff2(end:-1:1),theta);

    % Calculate the open-circuit potential of the active material
    OCP = refOCP + (T - Tref) .* dUdT;
    
end


%{
Copyright 2021-2023 SINTEF Industry, Sustainable Energy Technology
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
