function rho = pressure2density(p, T, MW)
%PRESSURE2DENSITY Converts gas pressure to mass density.
%   Assuming ideal gas behavior

con = physicalConstants();

rho = p .* MW ./(con.R * T);

end

