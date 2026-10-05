


load_experimental_data();

% --- 5. Plot the Nyquist diagram ---
figure('Name', 'EIS', 'Color', 'w');

% Plot -Im(Z) against Re(Z)
plot(Z_real, Z_imag, 'o', ...
    'LineWidth', 1.5, ...
    'MarkerSize', 2, ...
    'MarkerFaceColor', [0 0.4470 0.7410], ... 
    'MarkerEdgeColor', 'k');
grid on;
axis equal; % Use the same scale for the X and Y axes

% Add labels
xlabel('Z_{re} ', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('-Z_{im}', 'FontSize', 12, 'FontWeight', 'bold');
title('Nyquist diagram from data', 'FontSize', 14);

% Improve axis formatting
set(gca, 'FontSize', 11, 'LineWidth', 1);

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
