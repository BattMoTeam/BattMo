function c = molal2molar(b, M, rho)
%MOLAL2MOLAR Converts molality to molarity
%   Detailed explanation goes here

c = (b .* rho) ./ (1 + b .* M);

end

