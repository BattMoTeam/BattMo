function conductivity = computeElectrolyteConductivity_Xu2015(c, T)
    
    conductivityFactor = 1e-4;
    
    %cnst = [-10.5   , 0.074    , -6.96e-5; ...
    %        0.668e-3, -1.78e-5 , 2.80e-8; ...
    %        0.494e-6, -8.86e-10, 0];            
            
    
    % Ionic conductivity, [S m^-1]
    %conductivity = conductivityFactor.* c .*( polyval(cnst(end:-1:1,1),c) + polyval(cnst(end:-1:1,2),c) .* T + ...
    %                                          polyval(cnst(end:-1:1,3),c) .* T.^2).^2;
    %From cideMOD
    conductivity=c.*1e-4*1.2544.* ...
        (-8.2488+0.053248.*T-2.987e-5*(T.^2)+ ...
        0.26235e-3.*c-9.3063e-6.*c.*T+ ...
        8.069e-9.*c.*T.^2+ ...
        2.2002e-7.*c.^2-1.765e-10.*T.*c.^2);
    
end


%{
%}
