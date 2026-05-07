function [Z_real, Z_imag, omega] = load_experimental_data(filename)
    
    % --- 1. Configuration et Importation des données ---
    % filename = 'C:\Users\Alexandre Fichter\Documents\stage_3A\contenu stage\data_August\ank_data\Supplementary material\02_Electrical_characterization\EIS\131-828_EIS_01_MB_CD8.txt';
    
    % Le délimiteur est la tabulation ('\t')
    opts = detectImportOptions(filename, 'Delimiter', '\t');

    opts.VariableNamesLine = 1;
    opts.DataLine = 2;
    opts.VariableNamingRule = 'preserve';
    data = readtable(filename, opts);
    colNames = data.Properties.VariableNames;
    

    idx_real = contains(colNames, 'Re(Z)', 'IgnoreCase', true);
    idx_imag = contains(colNames, 'Im(Z)', 'IgnoreCase', true);
    idx_freq = contains(colNames, 'freq', 'IgnoreCase', true);
    
    % Extraction des données brutes (peu importe où elles se trouvent)
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
    Z_imag = Z_imag(idx_valides);
    freq = freq(idx_valides);

    Z_real = Z_real(8:60);
    Z_imag = Z_imag(8:60);
    omega = 2*pi*freq(8:60);

    Z_real = Z_real(:);
    Z_imag = Z_imag(:);
    omega = omega(:);
%% Figure test
    % figure;
    % semilogx(omega, Z_real, 'ro', 'MarkerFaceColor', 'r');      
    % % legend('Expérimental', 'Modèle (Fitted)');
    % title('Fitting results');
    % xlabel('Omega');
    % ylabel('Z_{re} '); 
    
end
