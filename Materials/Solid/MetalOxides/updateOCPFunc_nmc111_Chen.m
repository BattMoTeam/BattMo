function OCP = updateOCPFunc_nmc111_Chen(c, T, cmax)

    % LG M50 NMC open circuit potential as a function of stochiometry, fit taken
    % from [1].

    % [1] Chang-Hui Chen, Ferran Brosa Planella, Kieran O’Regan, Dominika Gastol, W.
    % Dhammika Widanage, and Emma Kendrick. "Development of Experimental Techniques for
    % Parameterization of Multi-scale Lithium-ion Battery Models." Journal of the
    % Electrochemical Society 167 (2020): 080534.
        
    OCP = -0.8090*c + 4.4875 - 0.0428*tanh(18.5138*(c - 0.5542)) - 17.7326*tanh(15.7890*(c - 0.3117)) + 17.5842*tanh(15.9308*(c - 0.3120));
    
end