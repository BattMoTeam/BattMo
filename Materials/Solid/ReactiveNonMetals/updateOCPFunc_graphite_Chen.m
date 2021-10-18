function OCP = updateOCPFunc_graphite_Chen(c, T, cmax)
    % LG M50 graphite open circuit potential as a function of stochiometry, fit taken
    % from [1].

    % References
    % ----------
    % .. [1] Chang-Hui Chen, Ferran Brosa Planella, Kieran O’Regan, Dominika Gastol, W.
    % Dhammika Widanage, and Emma Kendrick. "Development of Experimental Techniques for
    % Parameterization of Multi-scale Lithium-ion Battery Models." Journal of the
    % Electrochemical Society 167 (2020): 080534.

    OCP = 1.9793*exp(-39.3631*c) + 0.2482 - 0.0909*tanh(29.8538*(c - 0.1234)) - 0.04478*tanh(14.9159*(c - 0.2769)) - ...
          0.0205*tanh(30.4444*(c - 0.6103));
    
end