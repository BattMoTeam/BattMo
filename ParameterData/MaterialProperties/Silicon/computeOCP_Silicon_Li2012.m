function [OCP, dUdT] = computeOCP_Silicon_Li2012(c, T, cmax)
% From
% Juchuan Li, Xingcheng Xiao, Fuqian Yang, Mark W. Verbrugge, Yang-Tse Cheng
% Potentiostatic Intermittent Titration Technique for Electrodes Governed by Diffusion and Interfacial Reaction
% doi : 10.1021/jp207919q
    % if isa(c, 'ADI')
    %     c.val = real(c.val);
    %     fprintf('ADI - imag max: %e\n', max(abs(imag(value(c)))));
    % end
    stoc = c/cmax;
    stoc_vals = [0.00476555
                 0.00485132
                 0.00773118
                 0.01576
                 0.0207659
                 0.0379671
                 0.0622965
                 0.0907105
                 0.132342
                 0.193282
                 0.246094
                 0.297905
                 0.363947
                 0.451333
                 0.52145
                 0.605778
                 0.685034
                 0.744981
                 0.859803
                 0.943123
                ];

    OCP_vals = [1.37339
                1.37339
                1.20969
                1.01016
                0.861795
                0.713367
                0.593035
                0.511045
                0.441773
                0.369843
                0.303071
                0.264437
                0.233401
                0.207368
                0.196771
                0.152851
                0.126861
                0.100972
                0.0671234
                0.0385545
               ];
    OCP = interpTable(stoc_vals, OCP_vals, stoc);
    dUdT = 0;
    
end
