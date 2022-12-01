function res = butlerVolmer_split(j0, PF1, PF2, alpha, n, eta, T)
%BUTLERVOLMER Summary of this function goes here
%   Detailed explanation goes here
con = physicalConstants();

res = j0 .* (   PF1 .* exp(  alpha .* n .* con.F .* eta ./ ( con.R .* T ) ) - ...
                PF2 .* exp( -(1-alpha) .* n .* con.F .* eta ./ ( con.R .* T ) ) );


end

