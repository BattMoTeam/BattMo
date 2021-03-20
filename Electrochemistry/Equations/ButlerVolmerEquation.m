function res = ButlerVolmerEquation(j0, alpha, n, eta, T)
%BUTLERVOLMER Implements the standard form of the Butler-Volmer equation
%for electrochemical charge-transfer reaction kinetics.
%   Detailed explanation goes here
constants = PhysicalConstants();

res = j0 .* (   exp(  alpha .* n .* constants.F .* eta ./ ( constants.R .* T ) ) - ...
                exp( -(1-alpha) .* n .* constants.F .* eta ./ ( constants.R .* T ) ) );

end

