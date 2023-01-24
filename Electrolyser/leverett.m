function J = leverett(par, s)
    %LEVERETT calculates the leverett J-function of a porous
    %domain. 
    %   obj.leverett() calculates the leverett J-function
    %   for a liquid saturation of s, considering the properties of
    %   the porous medium.
    
    J0  = par(1);
    A1  = par(2);
    A2  = par(3);
    B1  = par(4);
    B2  = par(5);

    J   =   J0 + ...
            A1 .* exp(B1 .* (s-0.5)) - ...
            A2 .* exp(B2 .* (s-0.5));
                      
end
