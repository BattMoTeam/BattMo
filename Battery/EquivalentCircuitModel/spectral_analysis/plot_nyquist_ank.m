


load_experimental_data();

% --- 5. Tracé du diagramme de Nyquist ---
figure('Name', 'EIS', 'Color', 'w');

% On trace -Im(Z) en fonction de Re(Z)
plot(Z_real, Z_imag, 'o', ...
    'LineWidth', 1.5, ...
    'MarkerSize', 2, ...
    'MarkerFaceColor', [0 0.4470 0.7410], ... 
    'MarkerEdgeColor', 'k');
grid on;
axis equal; %  Force la même échelle en X et Y

% Ajout des labels
xlabel('Z_{re} ', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('-Z_{im}', 'FontSize', 12, 'FontWeight', 'bold');
title('Nyquist diagram from data', 'FontSize', 14);

% Amélioration des axes
set(gca, 'FontSize', 11, 'LineWidth', 1);