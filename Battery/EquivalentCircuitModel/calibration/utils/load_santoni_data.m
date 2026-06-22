function [Z_real, Z_imag, omegas] = load_santoni_data()

    path_data = 'C:\Users\Alexandre Fichter\Documents\stage_3A\contenu stage\papers\LiPo Battery LP-503562-IS-3 EIS, Capacity, ECM Data\LiPO_1\EIS_Charge_discharge\EIS_45\1_EIS.csv';
    data_exp = readtable(path_data, 'Delimiter', '\t');
    
   
    freq = table2array(data_exp(:, 1));
    Z_real = table2array(data_exp(:, 2));
    Z_imag = table2array(data_exp(:, 3));
    
    omegas = 2 * pi .* freq;
end