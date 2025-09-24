function OCP = computeOCP_Graphite_Latz(theta)
% computeOCP_Graphite_Latz
%   Computes the equilibrium open-circuit potential (OCP) of
%   graphite based on the model from Hein, Danner, and Latz [1].
%
%   Inputs:
%       theta - state of charge
%    
%   Output:
%       OCP   - open-circuit potential at 298.15 K [V vs Li/Li+]
%
%   Reference:
%       Hein, Danner, and Latz, ACS Appl. Energy Mater. 2020, 3, 8519−8531

    % Modified hyperbolic tangent function used in the original paper
    tanhmod = @(x) (exp(20 .* x) - exp(-x)) ./ (exp(-x) + exp(x));

    % OCP expression at reference temperature (T = 298.15 K)
    OCP = 0.6379 ...
        + 0.5416 .* exp(-305.5309 .* theta) ...
        + 0.044  .* tanh((-theta - 0.1958) ./ 0.1088) ...
        - 0.1978 .* tanhmod((theta - 0.99) ./ 0.05) ...
        - 0.6875 .* tanh((theta + 0.0117) ./ 0.0529) ...
        - 0.0175 .* tanh((theta - 0.5692) ./ 0.0875);


end
