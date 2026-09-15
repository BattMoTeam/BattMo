function R_lith = computeRadius(cAverage, cmax, R_delith)


            molarVolumeSi = 1.2e-05;
            molarVolumeLi = 8.8e-06;

            % We cannot anymore express <x> ( = theta) as a ratio of between the current and maximal concentrations of Li because the radius of the 
            % particle is changing. We have to express it as a ratio of matter quantities N/Nmax. The expression above is from a
            % rearrangement of equations 11 and 14 (can seem complex but can be done manually easily) in 'Analysis of Lithium Insertion/Deinsertion in a Silicon Electrode
            % Particle at Room Temperature Rajeswari, Chandrasekaran,Alexandre Magasinski, Gleb Yushin, and Thomas F. Fuller'*


            % defining useful quantities
            Q = (3.75.*molarVolumeLi)./(molarVolumeSi);
            theta = (cAverage./cmax);

            radius = R_delith.*(1 + Q.*theta).^(1/3);



            %Q = 3.75.* molarVolumeLi /(molarVolumeSi+3.75.*molarVolumeLi);
            %num = 1-Q;
            %denom = 1 - Q.* c_ratio;
%
            %radius = R_delith * ((3.85 .* num ./ denom).^ (1/3));


            if radius < R_delith
                error('R_lith has became smaller than R_delith; impossible')
            end

            R_lith = radius;

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
