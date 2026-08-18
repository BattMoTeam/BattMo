function conductivity = computeElectrolyteConductivity_Ai2020(c, T)
%    Conductivity of LiPF6 in EC:DMC as a function of ion concentration.
%
%    References
%    ----------
%    .. [1] Ai, W., Kraft, L., Sturm, J., Jossen, A., & Wu, B. (2020).
%    Electrochemical Thermal-Mechanical Modelling of Stress Inhomogeneity
%    in Lithium-Ion Pouch Cells. Journal of The Electrochemical Society,
%    167(1), 013512. DOI: 10.1149/2.0122001JES.
%    .. [2] Torchio, Marcello, et al. "Lionsimba: a matlab framework based
%    on a finite volume model suitable for li-ion battery design, simulation,
%    and control." Journal of The Electrochemical Society 163.7 (2016): A1192.
    
    conductivityFactor = 1e-4;
    
    cnst = [-10.5   , 0.074    , -6.96e-5; ...
            0.668e-3, -1.78e-5 , 2.80e-8; ...
            0.494e-6, -8.86e-10, 0];            
            
    
    % Ionic conductivity, [S m^-1]
    conductivity = conductivityFactor.* c .*( polyval(cnst(end:-1:1,1),c) + polyval(cnst(end:-1:1,2),c) .* T + ...
                                              polyval(cnst(end:-1:1,3),c) .* T.^2).^2;
    
end


%{
%}
