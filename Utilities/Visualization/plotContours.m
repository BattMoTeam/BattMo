function plotContours(model, states)

    assert(model.grid.griddim == 1, 'This function only works for 1D grids');

    time = cellfun(@(s) s.time, states) / hour;
    cmax_ne = model.NegativeElectrode.Coating.ActiveMaterial.Interface.saturationConcentration;
    cmax_pe = model.PositiveElectrode.Coating.ActiveMaterial.Interface.saturationConcentration;

    x_ne = model.NegativeElectrode.Coating.grid.cells.centroids/(micro*meter);
    x_pe = model.PositiveElectrode.Coating.grid.cells.centroids/(micro*meter);
    x_elyte = model.Electrolyte.grid.cells.centroids/(micro*meter);

    % Get the centroid locations and the values for the electrolyte
    % concentration at the given timestep
    c_ne = zeros(numel(states), numel(x_ne));
    c_pe = zeros(numel(states), numel(x_pe));
    c_elyte = zeros(numel(states), numel(x_elyte));
    phi_ne = zeros(numel(states), numel(x_ne));
    phi_pe = zeros(numel(states), numel(x_pe));
    phi_elyte = zeros(numel(states), numel(x_elyte));

    for istate = 1:numel(states)
        c_ne(istate, :) = states{istate}.NegativeElectrode.Coating.ActiveMaterial.SolidDiffusion.cSurface;
        c_pe(istate, :) = states{istate}.PositiveElectrode.Coating.ActiveMaterial.SolidDiffusion.cSurface;
        c_elyte(istate, :) = states{istate}.Electrolyte.c;
        phi_ne(istate, :) = states{istate}.NegativeElectrode.Coating.phi;
        phi_pe(istate, :) = states{istate}.PositiveElectrode.Coating.phi;
        phi_elyte(istate, :) = states{istate}.Electrolyte.phi;
    end

    % Plot the concentration values at the given grid centroid
    % locations
    figure
    contourf(x_ne, time, c_ne / cmax_ne, 20, 'LineStyle', 'none')
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
    contourf(x_pe, time, c_pe / cmax_pe, 20, 'LineStyle', 'none')
    xlabel('Position  /  µm')
    ylabel('Time  /  h')
    title('Positive Electrode Lithiation  /  -')
    colorbar()
    cm = crameri('nuuk');
    colormap(cm)
    set(gca, 'FontSize', 18);

    % Plot the potential at the given grid centroid locations
    figure
    contourf(x_ne, time, phi_ne, 20, 'LineStyle', 'none')
    xlabel('Position  /  µm')
    ylabel('Time  /  h')
    title('Negative Electrode Potential  /  V')
    colorbar()
    cm = flipud(crameri('lajolla'));
    colormap(cm)
    set(gca, 'FontSize', 18);

    figure
    contourf(x_elyte, time, phi_elyte, 20, 'LineStyle', 'none')
    xlabel('Position  /  µm')
    ylabel('Time  /  h')
    title('Electrolyte Potential  /  V')
    colorbar()
    cm = cmocean('curl'); %, 'pivot', 1);
    colormap(cm)
    set(gca, 'FontSize', 18);

    figure
    contourf(x_pe, time, phi_pe, 20, 'LineStyle', 'none')
    xlabel('Position  /  µm')
    ylabel('Time  /  h')
    title('Positive Electrode Potentials  /  V')
    colorbar()
    cm = crameri('nuuk');
    colormap(cm)
    set(gca, 'FontSize', 18);

end
