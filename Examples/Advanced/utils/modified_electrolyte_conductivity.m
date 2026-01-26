function conductivity = modified_electrolyte_conductivity(c, T)
%% modified version of electrolyte conductivity, used to enhance thermal coupling in an example
%% original is computeElectrolyteConductivity_Ai2020.m
    
    conductivityFactor = 1e-4;
    
    cnst = [-10.5   , 0.074    , -6.96e-5; ...
            0.668e-3, -1.78e-5 , 2.80e-8; ...
            0.494e-6, -8.86e-10, 0];            
            
    %
    modified_coef1 = 1; % original is 1.
    modified_coef2 = 1; % original is 1.
    
    % Ionic conductivity, [S m^-1]
    conductivity = conductivityFactor.* c .*( polyval(cnst(end:-1:1,1),c) + modified_coef1*polyval(cnst(end:-1:1,2),c) .* T + ...
                                                       modified_coef2*polyval(cnst(end:-1:1,3),c) .* T.^2).^2;
    % conductivity = max(conductivity, 1e-12);
    
end
