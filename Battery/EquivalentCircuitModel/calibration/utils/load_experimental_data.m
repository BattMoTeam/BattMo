function [Z_real, Z_imag, omega] = load_experimental_data(filename)
    
    % --- 1. Configuration and data import ---
    % filename = 'C:\Users\Alexandre Fichter\Documents\stage_3A\contenu stage\data_August\ank_data\Supplementary material\02_Electrical_characterization\EIS\131-828_EIS_01_MB_CD8.txt';
    
    % The delimiter is a tab ('\t')
    opts = detectImportOptions(filename, 'Delimiter', '\t');

    opts.VariableNamesLine = 1;
    opts.DataLine = 2;
    opts.VariableNamingRule = 'preserve';
    data = readtable(filename, opts);
    colNames = data.Properties.VariableNames;
    

    idx_real = contains(colNames, 'Re(Z)', 'IgnoreCase', true);
    idx_imag = contains(colNames, 'Im(Z)', 'IgnoreCase', true);
    idx_freq = contains(colNames, 'freq', 'IgnoreCase', true);
    
    % Extract raw data regardless of column position
    Z_real_raw = data{:, idx_real}; 
    Z_imag_raw = data{:, idx_imag};
    freq_raw  = data{:, idx_freq};

    if iscell(Z_real_raw) || isstring(Z_real_raw)
        Z_real = str2double(strrep(Z_real_raw, ',', '.'));
        Z_imag = str2double(strrep(Z_imag_raw, ',', '.'));
        freq  = str2double(strrep(freq_raw, ',', '.'));
    else
        Z_real = Z_real_raw;
        Z_imag = Z_imag_raw;    
        freq  = freq_raw;
    end

    idx_valides = (Z_real ~= 0) & ~isnan(Z_real);
    Z_real = Z_real(idx_valides);
    Z_imag = -Z_imag(idx_valides);                   % in the file is given -Im
    freq = freq(idx_valides);

    Z_real = Z_real(8:60);
    Z_imag = Z_imag(8:60);
    omega = 2*pi*freq(8:60);

    Z_real = Z_real(:);
    Z_imag = Z_imag(:);
    omega = omega(:);
%% Test plot
    % figure;
    % semilogx(omega, Z_real, 'ro', 'MarkerFaceColor', 'r');      
    % % legend('Experimental', 'Model (Fitted)');
    % title('Fitting results');
    % xlabel('Omega');
    % ylabel('Z_{re} '); 
    
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
