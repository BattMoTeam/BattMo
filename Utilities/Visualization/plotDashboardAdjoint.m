function [fig, tlo] = plotDashboardAdjoint(model, lambdas, varargin)

    opt = struct('step', 1, ...
                 'pauseTime', 0.1, ...
                 'scaleAxis', true);
    opt = merge_options(opt, varargin{:});

    % Shorthands
    ne      = 'NegativeElectrode';
    pe      = 'PositiveElectrode';
    co      = 'Coating';
    am      = 'ActiveMaterial';
    cc      = 'CurrentCollector';
    elyte   = 'Electrolyte';
    thermal = 'ThermalModel';
    ctrl    = 'Control';
    sd      = 'SolidDiffusion';

    % Parse the output from the lambdas
    eqnames = model.equationNames;
    varnames =  model.equationVarNames;
    assert(numel(eqnames) == numel(lambdas{1}));
    grids = cell(numel(eqnames), 1);

    for ieq = 1:numel(eqnames)
        eqname = eqnames{ieq};
        varname = varnames{ieq};

        % Skip if control is the varname
        if any(strcmp(varname, ctrl))
            continue
        end

        % Skip if solid diffusion variable (for now)
        if any(strcmp(varname, sd))
            continue
        end

        % Get grid from the model (as deep as possible)
        for k = 1:numel(varname)
            prop = varname(1:end-k);
            submodel = getfield(model, prop{:});

            if isprop(submodel, 'G') && ~isempty(submodel.G)
                grids{ieq} = submodel.grid();
                break
            end

        end

    end

    % Count non-empty eqs
    neq = sum(~cellfun(@isempty, grids));

    % Plot
    fig = figure;
    tlo = tiledlayout(2, ceil(neq/2));

    % Filter vars for plotting
    idx = ~cellfun(@isempty, grids);
    grids = grids(idx);
    eqnames = eqnames(idx);

    if opt.step <= 0
        % Plot all time steps

        for k = 1:numel(lambdas)

            sgtitle(sprintf('Adjoint variables at step %d', k));
            data = lambdas{k}(idx);

            if k == 1
                tile = plotDashboardFcn(fig, grids, data, eqnames, 'scaleAxis', opt.scaleAxis);
            else
                plotDashboardFcn(fig, grids, data, eqnames, 'scaleAxis', opt.scaleAxis, 'tile', tile);
            end

            pause(opt.pauseTime);

        end

    else

        % Plot single step
        sgtitle(sprintf('Adjoint variables at step %d', opt.step));
        data = lambdas{opt.step}(idx);
        plotDashboardFcn(fig, grids, data, eqnames, 'scaleAxis', opt.scaleAxis);

    end

end


function tile = plotDashboardFcn(fig, grids, data, titles, varargin)

    opt = struct('scaleAxis', true,  ...
                 'tile', []);
    opt = merge_options(opt, varargin{:});

    figure(fig);

    for k = 1:numel(grids)

        if isempty(opt.tile)
            tile{k} = nexttile;
        else
            tile{k} = opt.tile{k};
            axes(tile{k});
        end

        plotCellData(grids{k}, data{k});
        grid on;

        if isempty(opt.tile)
            drawnow;

            if opt.scaleAxis
                scaleAxis(gca, 'x', 1/micro);
                xlabel(gca, 'Position  /  µm')
            else
                xlabel(gca, 'Position  /  m')
            end

            title(titles{k}, 'Interpreter', 'none');
        end
    end

    if ~isempty(opt.tile)
        drawnow;
    end

end
