%{
    This is a script to create a gif .
%}

% Define the output GIF file name
gif_filename = 'Examples/Gifs/pouch_cell_layers.gif';

% Create the figure once and pass it into plotBatteryGrid
figure_handle = figure('Color', 'none', 'Renderer', 'opengl');  % Create and reuse this figure with no background

% Capture each layer's plot and add them to the GIF
jsonstruct_material = parseBattmoJson('Examples/jsondatafiles/sample_input.json');
jsonfilename = 'Examples/JsonDataFiles/geometryMultiLayerPouch.json';
jsonstruct_geometry = parseBattmoJson(jsonfilename);

% Define a fun color scheme (keeping bright pink and orange)
funColors = [
    0.7529, 0.7529, 0.7529;  % Silver
    0.2, 0.5, 0.2;  % Dark Green (Forest Green)
    0.9, 0.5, 0.1;  % Bright Orange
    0.1, 0.3, 0.8;  % Dark Blue
    0.9, 0.5, 0.1;  % Bright Orange
];

num_layers = 5;  % Number of layers in the pouch cell
layer_height = jsonstruct_material.NegativeElectrode.Coating.thickness + jsonstruct_material.PositiveElectrode.Coating.thickness + jsonstruct_material.Separator.thickness + jsonstruct_material.NegativeElectrode.CurrentCollector.thickness + jsonstruct_material.PositiveElectrode.CurrentCollector.thickness;
layer_height*5

for i = 1:num_layers
    % Set the number of layers in the geometry
    jsonstruct_geometry.Geometry.nLayers = i;

    fjv = flattenJsonStruct(jsonstruct_geometry);
    
    % Load additional JSON files for control and simulation parameters
    jsonfilename = fullfile('Examples', 'jsondatafiles', 'cc_discharge_control.json');
    jsonstruct_control = parseBattmoJson(jsonfilename);
    jsonfilename = fullfile('Examples', 'jsondatafiles', 'simulation_parameters.json');
    jsonstruct_simparams = parseBattmoJson(jsonfilename);
    
    % Merge JSON structures
    jsonstruct = mergeJsonStructs({jsonstruct_geometry, ...
                                    jsonstruct_material, ...
                                    jsonstruct_control, ...
                                    jsonstruct_simparams}, 'warn', false);
    
    model = setupModelFromJson(jsonstruct);  % Setup the model

    % Set current figure
    set(0, 'CurrentFigure', figure_handle);  % Ensure we use the same figure
    cla;  % Clear the axes, not the figure
    
    % Use the wrapper function to plot the battery model with your custom colors
    plotBatteryGridWithColors(model, 'fig', figure_handle, 'setstyle', false, 'colors', funColors);
    
    % Remove the legend if it gets created
    legend off;

    % Remove the axes and background
    axis off;  % Hide the axes
    set(gca, 'color', 'none');  % Set the axes background to transparent
    set(gcf, 'color', 'none');  % Set the figure background to transparent

    zlim([0 17e-4]);
    %axis tight

    view(45, 40)  % Set the view angle

    
    % Remove axes and set transparent background
    set(gca, 'Color', 'none');  % Ensure axes background is transparent
    set(gcf, 'Color', 'none');  % Ensure figure background is transparent

    drawnow  % Update the figure
    
    % Capture the frame and append to the GIF
    frame = getframe(gcf);  % Capture the current figure
    im = frame2im(frame);  % Convert the frame to an image
    im = im(:, :, 1:3);  % Remove alpha channel if exists
    [imind, cm] = rgb2ind(im, 256);  % Convert the frame to indexed image
    
    % Write the frame to the GIF
    if i == 1
        % First frame: create the GIF with infinite loop
        imwrite(imind, cm, gif_filename, 'gif', 'Loopcount', inf, 'DelayTime', 1);
    else
        % Subsequent frames: append to the GIF
        imwrite(imind, cm, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', 1);
    end
end

% Wrapper function to plot the battery grid with custom colors
function plotBatteryGridWithColors(model, varargin)
    % Default options for the wrapper function
    opt = struct('colors', [], ...
                 'setstyle', true, ...
                 'fig', [], ...
                 'shortLegendText', false, ...
                 'legendLocation', 'sw', ...
                 'axisLabels', false);
    opt = merge_options(opt, varargin{:});

    ne    = 'NegativeElectrode';
    pe    = 'PositiveElectrode';
    cc    = 'CurrentCollector';
    sep   = 'Separator';
    co    = 'Coating';

    % Use default colors if not provided
    if isempty(opt.colors)
        colors = crameri('vik', 5);  % Default color scheme if not specified
    else
        colors = opt.colors;  % Use the provided color scheme
    end

    if isempty(opt.fig)
        fig = figure;
    else
        fig = figure(opt.fig);
    end
    hold on

    G = model.grid;

    if G.griddim == 1
        edgeparams = {'linewidth', 2};
        facecolorname = 'color';
    else
        edgeparams = {'edgealpha', 0.5, 'edgecolor', 0*[1, 1, 1]};
        facecolorname = 'facecolor';
    end

    if model.include_current_collectors
        plotGrid(model.(pe).(cc).grid, facecolorname, colors(5,:), edgeparams{:});
    end

    plotGrid(model.(pe).(co).grid, facecolorname, colors(4,:), edgeparams{:});
    plotGrid(model.(sep).grid,     facecolorname, colors(3,:), edgeparams{:});
    plotGrid(model.(ne).(co).grid, facecolorname, colors(2,:), edgeparams{:});

    if model.include_current_collectors
        plotGrid(model.(ne).(cc).grid, facecolorname, colors(1,:), edgeparams{:});
    end

    if G.griddim == 3
        view(3)
    end

    legend off;  % Remove legend to avoid duplication
    axis tight;

    if opt.setstyle
        setFigureStyle('quantity', 'single');
    end

    drawnow();
end





