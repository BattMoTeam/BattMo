function figureHandle = plotGittEcmCalibration(result)
%PLOTGITTECMCALIBRATION Plot directional GITT calibration results.

    figureHandle = figure('Name', 'GITT ECM calibration');
    layout = tiledlayout(2, 2, 'TileSpacing', 'compact', ...
                         'Padding', 'compact');

    nexttile(layout);
    hold on
    plot(result.charge.ocvMap.soc, result.charge.ocvMap.voltage, ...
         'LineWidth', 1.5, 'DisplayName', 'Charge');
    plot(result.discharge.ocvMap.soc, result.discharge.ocvMap.voltage, ...
         'LineWidth', 1.5, 'DisplayName', 'Discharge');
    hold off
    xlabel('SOC / 1');
    ylabel('OCV / V');
    title('Directional OCV');
    legend('Location', 'best');
    grid on

    nexttile(layout);
    hold on
    plotDirectionalParameters(result, 'R0', '-');
    plotDirectionalParameters(result, 'R1', '--');
    plotDirectionalParameters(result, 'R2', ':');
    hold off
    set(gca, 'YScale', 'log');
    xlabel('SOC / 1');
    ylabel('Resistance / Ohm');
    title('2-RC resistances');
    legend('Location', 'best');
    grid on

    nexttile(layout);
    hold on
    plotDirectionalParameters(result, 'tau1', '-');
    plotDirectionalParameters(result, 'tau2', '--');
    hold off
    set(gca, 'YScale', 'log');
    xlabel('SOC / 1');
    ylabel('Time constant / s');
    title('2-RC time constants');
    legend('Location', 'best');
    grid on

    nexttile(layout);
    validation = result.validation;
    if isempty(validation)
        axis off
        text(0.5, 0.5, 'No held-out windows in full-calibration mode', ...
             'HorizontalAlignment', 'center');
    else
        hold on
        for direction = ["Charge", "Discharge"]
            mask = validation.Direction == direction;
            plot(validation.SOC(mask), 1e3*validation.RMSE2RC_V(mask), 'o-', ...
                 'LineWidth', 1.2, 'DisplayName', direction + " 2-RC");
            if any(isfinite(validation.RMSELegacy1RC_V(mask)))
                plot(validation.SOC(mask), ...
                     1e3*validation.RMSELegacy1RC_V(mask), 'x--', ...
                     'LineWidth', 1.0, 'DisplayName', direction + " legacy 1-RC");
            end
        end
        hold off
        xlabel('SOC / 1');
        ylabel('Time-weighted RMSE / mV');
        title('Held-out pulse windows');
        legend('Location', 'best');
        grid on
    end
end

function plotDirectionalParameters(result, field, lineStyle)
    directions = {'charge', 'discharge'};
    colors = lines(2);
    for i = 1:numel(directions)
        direction = directions{i};
        map = result.(direction).parameterMap;
        plot(map.SOC, map.(field), lineStyle, 'Color', colors(i, :), ...
             'LineWidth', 1.4, ...
             'DisplayName', sprintf('%s %s', capitalize(direction), field));
    end
end

function output = capitalize(input)
    output = [upper(input(1)), input(2:end)];
end
