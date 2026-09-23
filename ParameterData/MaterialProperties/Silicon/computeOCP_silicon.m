function [OCP, dUdT] = computeOCP_silicon(c, T, cmax)
% Calculate the equilibrium open-circuit potential of silicon according to the model
% given by Eq. (5) in Li et al. (2018) [1].
%
% References
% ----------
% .. [1] Li, H., Song, Y., Lu, B., and Zhang, J. (2018).
%    Effects of stress dependent electrochemical reaction on voltage hysteresis of
%    lithium ion batteries. Applied Mathematics and Mechanics, 39(10), 1453-1464.
%    DOI: 10.1007/s10483-018-2373-8.
% .. [2] Torchio, M., Magni, L., Gopaluni, R. B., Braatz, R. D., and Raimondo, D. M. (2016).
%    LIONSIMBA: A Matlab Framework Based on a Finite Volume Model Suitable for Li-Ion Battery
%    Design, Simulation, and Control. Journal of The Electrochemical Society, 163(7),
%    A1192-A1205. DOI: 10.1149/2.0291607jes.

% The optional thermal correction uses graphite coefficients from [2] as a proxy for silicon.

% It is a swelling Material --> theta cannot just be expressed by c/cmax,
% necessary to come back to the real definition in term of matter
% quantities : N/Nmax.


    use_graphite_thermal = false;


    c_ratio = c./cmax;

    R_delith      = 60e-09;

    molarVolumeSi = 1.2e-05;
    molarVolumeLi = 9e-06;

    Q = (3.75.*molarVolumeLi)./(molarVolumeSi);

    radius = computeRadius(c, cmax, R_delith);

    soc = c_ratio .* ((radius ./ R_delith).^3) ./(1 + Q);

    z = soc;

    % Calculate the open-circuit potential at the reference temperature for the given lithiation
    refOCP = (0.62 ...
              - 1.94 .* z ...
              + 5.8 .*  z.^2 ...
              - 7.13 .* z.^3 ...
              - 1.8 .*  z.^4 ...
              + 9.34 .* z.^5 ...
              - 4.76 .* z.^6);

    if use_graphite_thermal
        Tref = 298.15;  % [K]

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

        % Evaluate the graphite entropy fit at the swelling-corrected state of charge.
        dUdT = 1e-3 .* polyval(coeff1(end:-1:1), soc) ./ polyval(coeff2(end:-1:1), soc);
        OCP = refOCP + (T - Tref) .* dUdT;
    else
        dUdT = 0;
        OCP = refOCP;
    end

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
