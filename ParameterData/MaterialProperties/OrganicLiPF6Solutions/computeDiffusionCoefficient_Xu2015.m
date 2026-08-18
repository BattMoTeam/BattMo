function D = computeDiffusionCoefficient_Xu2015(c, T)
    
    % Calculate diffusion coefficients constant for the diffusion coefficient calculation
    cnst = [ -4.43, -54; 
             -0.22, 0.0 ];

    Tgi = [ 229; 5.0 ];
    
    % Diffusion coefficient, [m^2 s^-1]
    %Removed 10⁻⁴ otherwise the same
    D = 10 .^ ( ( cnst(1,1) + cnst(1,2) ./ ( T - Tgi(1) - Tgi(2) .* c .* 1e-3) + cnst(2,1) .* ...
                          c .* 1e-3) );
    
end



%{
%}
