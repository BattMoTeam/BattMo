function res = butlerVolmer(j0, alpha, n, eta, T)
%BUTLERVOLMER Summary of this function goes here
%   Detailed explanation goes here
con = physicalConstants();

res = j0 .* (   exp(  alpha .* n .* con.F .* eta ./ ( con.R .* T ) ) - ...
                exp( -(1-alpha) .* n .* con.F .* eta ./ ( con.R .* T ) ) );


end

