
%  params = [6.2261e-3, 0.3e-3, 10000.1264, 0.00097, 2e5];  % best params
% handmade

params = [0.003, 0.00134, 5, 0.00330, 20000];
[Z_re_exp, Z_im_exp, omega] = load_experimental_data();
omega = logspace(-2, 4, 100);

[Z_real, Z_imag] = load_nyquist(params, omega)

figure;
plot(Z_re_exp, Z_im_exp, 'ro', 'MarkerFaceColor', 'r');
hold on;
plot(Z_real, Z_imag, '-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
axis equal; % Essential for a Nyquist plot
grid on;
xlabel('Z_{réel} (\Omega)', 'FontWeight', 'bold');
ylabel('-Z_{imaginaire} (\Omega)', 'FontWeight', 'bold');
title('Nyquist diagram simulated');

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
