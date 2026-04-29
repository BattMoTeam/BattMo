function [Z_real, Z_imag, omega] = load_experimental_data()
    % --- 1. Configuration et Importation des données ---
    filename = 'C:\Users\Alexandre Fichter\Documents\stage_3A\contenu stage\data_August\ank_data\Supplementary material\02_Electrical_characterization\EIS\131-828_EIS_01_MB_CD8.txt';
    
    % Le délimiteur est la tabulation ('\t')
    opts = detectImportOptions(filename, 'Delimiter', '\t');
    data = readtable(filename, opts);
    
    % --- 2. Identification des colonnes ---
    % L'impédance réelle est dans Var26 et l'imaginaire dans Var27.
    Z_real_raw = data.Var26; 
    Z_imag_raw = data.Var27;
    omega_raw = data.Var16;
    % --- 3. Gestion des virgules ---
    % On remplace les virgules par des points, puis on convertit le texte en nombres
    Z_real = str2double(strrep(Z_real_raw, ',', '.'));
    Z_imag = str2double(strrep(Z_imag_raw, ',', '.'));
    omega = str2double(strrep(omega_raw, ',', '.'));

    % --- 4. Nettoyage des données ---
    % S'il y a des lignes de cyclage sans données EIS, on filtre les zéros ou NaN
    idx_valides = (Z_real ~= 0) & ~isnan(Z_real);
    Z_real = Z_real(idx_valides);
    Z_imag = Z_imag(idx_valides);
    omega = omega(idx_valides);

    Z_real = Z_real(8:60);
    Z_imag = Z_imag(8:60);
    omega = 2*pi*omega(8:60);

    Z_real = Z_real(:);
    Z_imag = Z_imag(:);
    omega = omega(:);
    
end