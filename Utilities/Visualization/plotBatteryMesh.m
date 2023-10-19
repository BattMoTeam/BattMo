function plotBatteryMesh(model, setstyle, fig)

    if nargin == 1
        setstyle = true;
        fig = [];
    elseif nargin == 2
        fig = [];
    end

    ne    = 'NegativeElectrode';
    pe    = 'PositiveElectrode';
    cc    = 'CurrentCollector';
    sep   = 'Separator';
    co    = 'Coating';

    colors = crameri('vik', 5);
    if isempty(fig)
        figure
    else
        figure(fig);
    end
    hold on
    legtext = {};

    if model.G.griddim == 1
        edgeparams = {'linewidth', 2};
        facecolorname = 'color';
    else
        edgeparams = {'edgealpha', 0.5, 'edgecolor', 0*[1, 1, 1]};
        facecolorname = 'facecolor';
    end

    if model.include_current_collectors
        plotGrid(model.(pe).(cc).G, facecolorname, colors(5,:), edgeparams{:});
        legtext{end+1} = 'positive electrode current collector';
    end

    plotGrid(model.(pe).(co).G, facecolorname, colors(4,:), edgeparams{:});
    plotGrid(model.(sep).G,     facecolorname, colors(3,:), edgeparams{:});
    plotGrid(model.(ne).(co).G, facecolorname, colors(2,:), edgeparams{:});
    legtext = [legtext, {'positive electrode active material', 'separator', 'negative electrode active material'}];

    if model.include_current_collectors
        plotGrid(model.(ne).(cc).G, facecolorname, colors(1,:), edgeparams{:});
        legtext{end+1} = 'negative electrode current collector';
    end

    if model.G.griddim == 3
        view(3)
    end

    legend(legtext, 'location', 'southwest');
    axis tight;

    if setstyle
        setFigureStyle('quantity', 'single');
    end

    drawnow();

end
