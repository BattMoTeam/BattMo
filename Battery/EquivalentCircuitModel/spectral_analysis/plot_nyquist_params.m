
%  params = [6.2261e-3, 0.3e-3, 10000.1264, 0.00097, 2e5];  % best params
% handmade

params = [0.003, 0.00134, 5, 0.00330, 20000]
[Z_re_exp, Z_im_exp, omega] = load_experimental_data();
omega = logspace(-2, 4, 100)

[Z_real, Z_imag] = load_nyquist(params, omega)

figure;
plot(Z_re_exp, Z_im_exp, 'ro', 'MarkerFaceColor', 'r');
hold on;
plot(Z_real, Z_imag, '-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
axis equal; % Essentiel pour Nyquist
grid on;
xlabel('Z_{réel} (\Omega)', 'FontWeight', 'bold');
ylabel('-Z_{imaginaire} (\Omega)', 'FontWeight', 'bold');
title('Nyquist diagram simulated');







