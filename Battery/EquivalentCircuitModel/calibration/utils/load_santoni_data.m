function [Z_real, Z_imag, omegas] = load_santoni_data()

    path_data = 'C:\Users\Alexandre Fichter\Documents\stage_3A\contenu stage\papers\LiPo Battery LP-503562-IS-3 EIS, Capacity, ECM Data\LiPO_1\EIS_Charge_discharge\EIS_45\1_EIS.csv';
    data_exp = readtable(path_data, 'Delimiter', '\t');
    
   
    freq = table2array(data_exp(:, 1));
    Z_real = table2array(data_exp(:, 2));
    Z_imag = table2array(data_exp(:, 3));
    
    omegas = 2 * pi .* freq;
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
