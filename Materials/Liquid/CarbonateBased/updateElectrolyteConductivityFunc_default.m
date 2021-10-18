function conductivity = updateElectrolyteConductivityFunc_default(c, T)
    
    conductivityFactor = 1e-4;
    
    cnst = [-10.5   , 0.074    , -6.96e-5; ...
            0.668e-3, -1.78e-5 , 2.80e-8; ...
            0.494e-6, -8.86e-10, 0];            
            
    
    % Ionic conductivity, [S m^-1]
    conductivity = conductivityFactor.* c .*( polyval(cnst(end:-1:1,1),c) + polyval(cnst(end:-1:1,2),c) .* T + ...
                                              polyval(cnst(end:-1:1,3),c) .* T.^2).^2;
    
end