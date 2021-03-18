function res = ButlerVolmerEquation(j0, alpha, n, eta, T)
%BUTLERVOLMER Summary of this function goes here
%   Detailed explanation goes here
constants = PhysicalConstants();

res = j0 .* (   exp(  alpha .* n .* constants.F .* eta ./ ( constants.R .* T ) ) - ...
                exp( -(1-alpha) .* n .* constants.F .* eta ./ ( constants.R .* T ) ) );


end

