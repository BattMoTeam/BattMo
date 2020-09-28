function c = pressure2concentration(p, T)
%PRESSURE2CONCENTRATION Converts gas pressure to molar concentration
%   Assuming ideal gas behavior

con = physicalConstants();

c = p ./ (con.R .* T);

end
