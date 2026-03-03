function plotContours(output)

    time = output.time/hour;
    cmax_ne = output.model.NegativeElectrode.Coating.ActiveMaterial.Interface.saturationConcentration;
    cmax_pe = output.model.PositiveElectrode.Coating.ActiveMaterial.Interface.saturationConcentration;

    x_ne = output.model.NegativeElectrode.grid.cells.centroids/(micro*meter);
    x_pe = output.model.PositiveElectrode.grid.cells.centroids/(micro*meter);
    x_elyte = output.model.Electrolyte.grid.cells.centroids/(micro*meter);

    % Get the centroid locations and the values for the electrolyte
    % concentration at the given timestep

    c_elyte = cellfun(@(s) s.Electrolyte.c, output.states, 'UniformOutput', false);
    c_ne = cellfun(@(s) s.NegativeElectrode.Coating.ActiveMaterial.SolidDiffusion.cSurface, output.states, 'UniformOutput', false);
    c_pe = cellfun(@(s) s.PositiveElectrode.Coating.ActiveMaterial.SolidDiffusion.cSurface, output.states, 'UniformOutput', false);
    phi_elyte = cellfun(@(s) s.Electrolyte.phi, output.states, 'UniformOutput', false);
    phi_ne = cellfun(@(s) s.NegativeElectrode.Coating.phi, output.states, 'UniformOutput', false);
    phi_pe = cellfun(@(s) s.PositiveElectrode.Coating.phi, output.states, 'UniformOutput', false);

    % Convert cell arrays to matrices if needed
    c_elyte = cell2mat(c_elyte')';
    c_ne = cell2mat(c_ne')';
    c_pe = cell2mat(c_pe')';
    phi_elyte = cell2mat(phi_elyte')';
    phi_ne = cell2mat(phi_ne')';
    phi_pe = cell2mat(phi_pe')';

    % Plot the concentration values at the given grid centroid
    % locations
    figure
    contourf(x_ne, time, c_ne ./ cmax_ne, 20, 'LineStyle', 'none')
    xlabel('Position  /  µm')
    ylabel('Time  /  h')
    title('Negative Electrode Lithiation  /  -')
    colorbar()
    cm = flipud(crameri('lajolla'));
    colormap(cm)
    set(gca, 'FontSize', 18);

    figure
    contourf(x_elyte, time, c_elyte, 20, 'LineStyle', 'none')
    xlabel('Position  /  µm')
    ylabel('Time  /  h')
    title('Electrolyte Concentration  /  mol \cdot m^{-3}')
    colorbar()
    cm = cmocean('curl', 'pivot', 1000);
    colormap(cm)
    set(gca, 'FontSize', 18);

    figure
    contourf(x_pe, time, c_pe ./ cmax_pe, 20, 'LineStyle', 'none')
    xlabel('Position  /  µm')
    ylabel('Time  /  h')
    title('Positive Electrode Lithiation  /  -')
    colorbar()
    cm = crameri('nuuk');
    colormap(cm)
    set(gca, 'FontSize', 18);

    % % Plot the potential at the given grid centroid locations
    % figure
    % contourf(x_ne, time, phi_ne, 20, 'LineStyle', 'none')
    % xlabel('Position  /  µm')
    % ylabel('Time  /  h')
    % title('Negative Electrode Potential  /  V')
    % colorbar()
    % cm = flipud(crameri('lajolla'));
    % colormap(cm)
    % set(gca, 'FontSize', 18);

    % figure
    % contourf(x_elyte, time, phi_elyte, 20, 'LineStyle', 'none')
    % xlabel('Position  /  µm')
    % ylabel('Time  /  h')
    % title('Electrolyte Potential  /  V')
    % colorbar()
    % cm = cmocean('curl'); %, 'pivot', 1);
    % colormap(cm)
    % set(gca, 'FontSize', 18);

    % figure
    % contourf(x_pe, time, phi_pe, 20, 'LineStyle', 'none')
    % xlabel('Position  /  µm')
    % ylabel('Time  /  h')
    % title('Positive Electrode Potentials  /  V')
    % colorbar()
    % cm = crameri('nuuk');
    % colormap(cm)
    % set(gca, 'FontSize', 18);

end
