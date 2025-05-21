function [OCP] = computeOCP_Graphite_Latz(c, cmax)
% computeOCP_Graphite_Latz
%   Computes the equilibrium open-circuit potential (OCP) of
%   graphite based on the model from Hein, Danner, and Latz [1].
%
%   Inputs:
%       c     - lithium concentration in the solid phase [mol/m^3]
%       cmax  - maximum lithium concentration in the solid [mol/m^3]
%
%   Output:
%       OCP   - open-circuit potential at 298.15 K [V vs Li/Li+]
%
%   Reference:
%       Hein, Danner, and Latz, ACS Appl. Energy Mater. 2020, 3, 8519−8531

    % state of charge
    x = c ./ cmax;

    % Modified hyperbolic tangent function used in the original paper
    tanhmod = @(x) (exp(20 .* x) - exp(-x)) ./ (2 .* cosh(x));

    % OCP expression at reference temperature (T = 298.15 K)
    OCP = 0.6379 ...
        + 0.5416 .* exp(-305.5309 .* x) ...
        + 0.044  .* tanh((-x - 0.1958) ./ 0.1088) ...
        - 0.1978 .* tanhmod((x - 0.99) ./ 0.05) ...
        - 0.6875 .* tanh((x + 0.0117) ./ 0.0529) ...
        - 0.0175 .* tanh((x - 0.5692) ./ 0.0875);
end