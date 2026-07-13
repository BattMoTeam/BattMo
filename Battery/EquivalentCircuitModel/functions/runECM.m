function output = runECM(params, input, options)

    
    R0 = params(1);
    R1 = params(2);
    C1 = params(3);
    R2 = params(4);
    C2 = params(5);


    N_points = length(input.time_vec);
    V_ecm = zeros(N_points, 1);

    Uc1 = 0; 
    Uc2 = 0;
    
    % Sensibilities initialisation 
    
    s_R1 = 0; % dUc1 / dR1
    s_C1 = 0; % dUc1 / dC1
    s_R2 = 0; % dUc2 / dR2
    s_C2 = 0; % dUc2 / dC2
    
    g_true = zeros(5, 1);
    
    for k = 1:N_points
        
        if k == 1
            dt = input.time_vec(1);
        else
            dt = input.time_vec(k) - input.time_vec(k-1);
        end
        
        I = input.current_exp(k);           

        V_ocv = input.ocv_vec(k); 
        
        % 1. Tension du modèle au pas k
        V_ecm(k) = V_ocv - R0 * I - Uc1 - Uc2;
        
        % 2. Erreur instantanée
        err = input.voltage_exp(k) - V_ecm(k);

        if nargout > 1
            % Dérivées partielles de V_ecm par rapport à chaque paramètre
            dV_dR0 = -I;
            dV_dR1 = -s_R1;
            dV_dC1 = -s_C1;
            dV_dR2 = -s_R2;
            dV_dC2 = -s_C2;
            
            % Accumulation du gradient de la fonction coût (v = sum(err^2))
            g_true(1) = g_true(1) - 2 * err * dV_dR0;
            g_true(2) = g_true(2) - 2 * err * dV_dR1;
            g_true(3) = g_true(3) - 2 * err * dV_dC1;
            g_true(4) = g_true(4) - 2 * err * dV_dR2;
            g_true(5) = g_true(5) - 2 * err * dV_dC2;
            
            % Mise à jour des sensibilités pour le pas suivant (Euler Explicite)
            s_R1 = s_R1 + (-s_R1 / (R1 * C1) + Uc1 / (R1^2 * C1)) * dt;
            s_C1 = s_C1 + (-s_C1 / (R1 * C1) + Uc1 / (R1 * C1^2) - I / C1^2) * dt;
            s_R2 = s_R2 + (-s_R2 / (R2 * C2) + Uc2 / (R2^2 * C2)) * dt;
            s_C2 = s_C2 + (-s_C2 / (R2 * C2) + Uc2 / (R2 * C2^2) - I / C2^2) * dt;
        end
        
        % 4. Mise à jour des variables d'état électriques pour le pas suivant
        Uc1 = Uc1 + (-Uc1 / (R1 * C1) + I / C1) * dt;
        Uc2 = Uc2 + (-Uc2 / (R2 * C2) + I / C2) * dt;
    end
    
    % Valeur de la fonction objectif (Erreur quadratique totale)
    v = sum((input.voltage_exp - V_ecm).^2);
    

    if isnan(v) || isinf(v) || (nargout > 1 && (any(isnan(g_true)) || any(isinf(g_true))))
        v = 1e10; % Assigne un coût massif pour rejeter le point instable
        if nargout > 1
            g_norm = zeros(5, 1); % Renvoie un gradient plat pour forcer le pas arrière
        end
        return;
    end

    % 5. Application de la règle de dérivation en chaîne pour l'espace normalisé
    if nargout > 1
        pmin = input.scales(1:5);
        pmax = input.scales(6:10);
        dp_dpnorm = (pmax - pmin);      
        g_norm = g_true .* dp_dpnorm(:);
    end

end
