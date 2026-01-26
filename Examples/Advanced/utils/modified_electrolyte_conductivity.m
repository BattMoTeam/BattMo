function conductivity = modified_electrolyte_conductivity(c, T)
%% modified version of electrolyte conductivity, used to enhance thermal coupling in an example
    
    conductivityFactor = 1;
    
    cnst = [-10.5   , 0.074    , -6.96e-5; ...
            0.668e-3, -1.78e-5 , 2.80e-8; ...
            0.494e-6, -8.86e-10, 0];            
            
    
    % Ionic conductivity, [S m^-1]
    conductivity = conductivityFactor.* c .*( polyval(cnst(end:-1:1,1),c) + 1e1*polyval(cnst(end:-1:1,2),c) .* T + ...
                                              1e1*polyval(cnst(end:-1:1,3),c) .* T.^2).^2;
    
end
