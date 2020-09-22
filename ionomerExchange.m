function R = ionomerExchange(kxch, cT, z, inmrPhi, elytePhi, T, elyteC)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

con =   physicalConstants();

R   =   kxch .* ( cT .* ...
                ( exp( z.*con.F.*(inmrPhi - elytePhi) ./ (con.R .* T) ) )... -  ...
                  ...exp( -z.*con.F.*(inmrPhi - elytePhi) ./  (con.R .* T) ) ) ...
                - elyteC );
end

