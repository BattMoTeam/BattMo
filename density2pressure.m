function p = density2pressure(rho, T, MW)
%DENSITY2PRESSURE Converts gas pressure to mass density.
%   Assuming ideal gas behavior

con = physicalConstants();

p = rho .*(con.R * T) ./ MW ;

end