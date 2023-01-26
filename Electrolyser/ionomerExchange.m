function exchR = ionomerExchange(kxch, z, cInmr, cElyte, phiInmr, phiElyte, T)
%IONOMEREXCHANGE Calculates the ionomer to liquid ion-exchange rate
%assuming a first-order process

%   Follows an approach from Stanislaw, Gerhardt, and Weber ECS Trans 2019,
%   as well as Jiangjin Liu et al JES 2021.
%   The equation for R approaches the equation for Donnan equilibrium if
%   kxch is large or R approaches zero. So using a large kxch value should
%   enforce Donnan equilibrium.

%   kxch        Ionomer to liquid ion-exchange rate constant, 1/s. Best thing to do is
%               increase it until the converged solution no longer changes.
%   z           Ion charge
%   cInmr       Concentration of the relevant ion in the ionomer, mol/m^3. For a perfect
%               membrane with only one counter-ion, this should be equal to the
%               membrane fixed charge.
%   cElyte      Concentration of the relevant ion in the electrolyte,
%               mol/m^3.
%   Phiinmr     Ionic potential in the ionomer, V
%   Phielyte    Ionic potential in the electrolyte, V
%   T           Temperature, K

    con = PhysicalConstants();
    R = con.R;
    F = con.F;
    
    exchR = kxch.*(cInmr.*(exp(z*F*(phiInmr - phiElyte)./(R*T))) - cElyte);
    
end

