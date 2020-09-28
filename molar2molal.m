function b = molar2molal(c, M, rho)
%MOLAR2MOLAL Converts molarity to molality
%   Detailed explanation goes here

b = (c) ./ (rho - c .* M);

end

