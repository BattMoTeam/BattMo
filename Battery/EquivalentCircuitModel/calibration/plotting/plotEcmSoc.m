function plotEcmSoc(ecm_table)

    % Create a clean figure window
    figure('Color', 'w', 'Name', 'ECM Parameters vs SOC', 'Position', [100, 100, 1100, 700]);

    % 1. Plot Ohmic Resistance R0
    subplot(2,3,1);
    plot(ecm_table.SOC * 100, ecm_table.R0, 'b-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
    grid on;
    xlabel('SOC (%)'); ylabel('R_0 (\Omega)');
    title('Internal Resistance (R_0)');

    % 2. Plot RC1 Resistance R1
    subplot(2,3,2);
    plot(ecm_table.SOC * 100, ecm_table.R1, 'r-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'r');
    grid on;
    xlabel('SOC (%)'); ylabel('R_1 (\Omega)');
    title('Charge Transfer Resistance (R_1)');

    % 3. Plot RC1 Capacitance C1
    subplot(2,3,3);
    plot(ecm_table.SOC * 100, ecm_table.C1, 'm-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'm');
    grid on;
    xlabel('SOC (%)'); ylabel('C_1 (F)');
    title('Double Layer Capacitance (C_1)');

    % 4. Plot RC2 Resistance R2
    subplot(2,3,4);
    plot(ecm_table.SOC * 100, ecm_table.R2, 'g-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'g');
    grid on;
    xlabel('SOC (%)'); ylabel('R_2 (\Omega)');
    title('Diffusion Resistance (R_2)');

    % 5. Plot RC2 Capacitance C2
    subplot(2,3,5);
    plot(ecm_table.SOC * 100, ecm_table.C2, 'k-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'k');
    grid on;
    xlabel('SOC (%)'); ylabel('C_2 (F)');
    title('Diffusion Capacitance (C_2)');

    % Global title for the entire figure
    sgtitle('Fitted ECM Parameters Evolution across SOC Range', ...
            'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'none');
end
