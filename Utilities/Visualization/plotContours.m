function plotContours(model, states, varargin)

    assert(model.grid.griddim == 1, 'This function only works for 1D grids');

    opt = struct('subplot', true, ...
                 'varname', {'lithiation'}, ... % Or {'concentration'}
                 'sgtitle', '', ...
                 'plot_separators', false); % Not implemented
    opt = merge_options(opt, varargin{:});

    assert(containsi(opt.varname, 'lith') || containsi(opt.varname, 'conc'), ...
           'opt.varname must contain either ''lith'' or ''conc''');

    fontsize = 18;

    if opt.subplot
        fig = figure;
    end

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
    if opt.subplot
        figure(fig);
        subplot(2, 3, 1);
    else
        figure;
    end
    if containsi(opt.varname, 'lith')
        contourf(x_ne, time, c_ne / cmax_ne, 20, 'LineStyle', 'none')
        titlestr = 'Lithiation  /  -';
    elseif containsi(opt.varname, 'conc')
        contourf(x_ne, time, c_ne, 20, 'LineStyle', 'none')
        titlestr = 'Concentration  /  mol \cdot m^{-3}';
    end
    xlabel('Position  /  µm')
    ylabel('Time  /  h')
    title(sprintf('NE %s', titlestr));
    colorbar()
    cm = flipud(crameri('lajolla'));
    colormap(gca, cm)
    set(gca, 'FontSize', fontsize);
    drawnow

    if opt.subplot
        figure(fig);
        subplot(2, 3, 2);
    else
        figure;
    end
    contourf(x_elyte, time, c_elyte, 20, 'LineStyle', 'none')
    xlabel('Position  /  µm')
    ylabel('Time  /  h')
    title('Elyte Concentration  /  mol \cdot m^{-3}')
    colorbar()
    cm = cmocean('curl');
    cmscaled = rescaleCmap(cm, min(c_elyte), max(c_elyte), model.Electrolyte.species.nominalConcentration);
    colormap(gca, cmscaled)

    set(gca, 'FontSize', fontsize);
    drawnow

    if opt.subplot
        figure(fig);
        subplot(2, 3, 3);
    else
        figure;
    end
    if containsi(opt.varname, 'lith')
        contourf(x_pe, time, c_pe / cmax_pe, 20, 'LineStyle', 'none')
        titlestr = 'Lithiation  /  -';
    elseif containsi(opt.varname, 'conc')
        contourf(x_pe, time, c_pe, 20, 'LineStyle', 'none')
        titlestr = 'Concentration  /  mol \cdot m^{-3}';
    end
    xlabel('Position  /  µm')
    ylabel('Time  /  h')
    title(sprintf('PE %s', titlestr));
    colorbar()
    cm = crameri('nuuk');
    colormap(gca, cm)
    set(gca, 'FontSize', fontsize);
    drawnow

    % Plot the potential at the given grid centroid locations
    if opt.subplot
        figure(fig);
        subplot(2, 3, 4);
    else
        figure;
    end
    contourf(x_ne, time, phi_ne, 20, 'LineStyle', 'none')
    xlabel('Position  /  µm')
    ylabel('Time  /  h')
    title('NE Potential  /  V')
    colorbar()
    cm = flipud(crameri('lajolla'));
    colormap(gca, cm)
    set(gca, 'FontSize', fontsize);
    drawnow

    if opt.subplot
        figure(fig);
        subplot(2, 3, 5);
    else
        figure;
    end
    contourf(x_elyte, time, phi_elyte, 20, 'LineStyle', 'none')
    xlabel('Position  /  µm')
    ylabel('Time  /  h')
    title('Elyte Potential  /  V')
    colorbar()
    cm = cmocean('curl'); %, 'pivot', 1);
    colormap(gca, cm)
    set(gca, 'FontSize', fontsize);
    drawnow

    if opt.subplot
        figure(fig);
        subplot(2, 3, 6);
    else
        figure;
    end
    contourf(x_pe, time, phi_pe, 20, 'LineStyle', 'none')
    xlabel('Position  /  µm')
    ylabel('Time  /  h')
    title('PE Potentials  /  V')
    colorbar()
    cm = crameri('nuuk');
    colormap(gca, cm)
    set(gca, 'FontSize', fontsize);
    drawnow

    if opt.subplot
        if ~isempty(opt.sgtitle)
            sgtitle(opt.sgtitle, 'FontSize', fontsize+2);
        end
    end

end


function s = containsi(a, b)
    s = contains(a, b, 'IgnoreCase', true);
end


function cmap_scaled = rescaleCmap(cmap, vmin, vmax, center, Nout)

    if nargin < 5
        Nout = 256;
    end

    % Split colormap into two halves
    N = size(cmap,1);
    cmap_low  = cmap(1:floor(N/2),:);   % lower half (below center)
    cmap_high = cmap(floor(N/2)+1:end,:); % upper half (above center)

    % Number of colors for each side proportional to data range
    n_low  = round(Nout * (center - vmin) / (vmax - vmin));
    n_high = Nout - n_low;

    % Interpolate each side independently
    cmap_scaled = [ ...
        interp1(linspace(0,1,size(cmap_low,1)), cmap_low,  linspace(0,1,n_low)); ...
        interp1(linspace(0,1,size(cmap_high,1)), cmap_high, linspace(0,1,n_high)) ...
                  ];
end
