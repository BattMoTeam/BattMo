%% Fitting Script

% Charge les données expérimentales
t_exp = [0:1:1000]'; % Temps en secondes
V_exp = 4.1 * ones(1001, 1) - 0.1 * (t_exp/1000); % Fausse courbe de tension
I_exp = ones(1001, 1); % Courant expérimental


parameters = createParametersECM();

% Remplacer le courant du modèle par le courant expérimental
parameters.I.dataX = t_exp;
parameters.I.dataY = I_exp;

% Valeurs initiales [R0, R1, C1, R2, C2]
initial_values = [0.010, 0.010, 3000, 0.010, 100000];

options = optimset('Display', 'iter', 'TolX', 1e-4);
disp('beginning of optimization');

% Optimization with fminsearch
best_params = fminsearch(@(p) costFunctionECM(p, parameters, initial_values, t_exp, V_exp), initial_values, options);

fprintf('\n=== PARAMETERS FOUND ===\n');
fprintf('R0 = %.5f Ohms\n', best_params(1));
fprintf('R1 = %.5f Ohms\n', best_params(2));
fprintf('C1 = %.1f Farads\n', best_params(3));
fprintf('R2 = %.5f Ohms\n', best_params(4));
fprintf('C2 = %.1f Farads\n', best_params(5));

%% Vérification Graphique 


parameters.R0 = best_params(1);
parameters.R1 = best_params(2);
parameters.C1 = best_params(3);
parameters.R2 = best_params(4);
parameters.C2 = best_params(5);

inputparams = EquivalentCircuitModelInputParams(parameters);
model = EquivalentCircuitModel(inputparams);
[t_sim, V_sim, ~] = model.solve();

figure;
plot(t_exp, V_exp, 'k', 'LineWidth', 2); hold on;
plot(t_sim, V_sim, 'r--', 'LineWidth', 2);
legend('Expérimental', 'Modèle (Fitted)');
title('Résultat du Fitting');

%{
Copyright 2021-2026 SINTEF Industry, Sustainable Energy Technology
and SINTEF Digital, Mathematics & Cybernetics.

This file is part of The Battery Modeling Toolbox BattMo

BattMo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BattMo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BattMo.  If not, see <http://www.gnu.org/licenses/>.
%}
