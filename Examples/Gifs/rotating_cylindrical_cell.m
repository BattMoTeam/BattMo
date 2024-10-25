%{
    This is a script to create a gif .
%}

% Define the output GIF file name
gif_filename = 'Examples/Gifs/rotating_cylindrical_cell.gif';

% Number of rotation steps (higher value = smoother rotation)
num_frames = 36;  % This will give 360 degrees / 36 = 10 degrees per step

% Create the figure once and pass it into plotBatteryGrid
figure_handle = figure('Color', 'none');  % Create and reuse this figure with no background

% Load JSON files and setup the model
jsonstruct_material = parseBattmoJson('Examples/jsondatafiles/sample_input.json');
jsonstruct_material.NegativeElectrode.Coating.thickness = 1e-4
jsonstruct_material.PositiveElectrode.Coating.thickness = 1e-4
jsonfilename = 'Examples/JsonDataFiles/4680-geometry.json';
jsonstruct_geometry = parseBattmoJson(jsonfilename);
jsonstruct_geometry.Geometry.rOuter = jsonstruct_geometry.Geometry.rInner + 4*milli*meter;
jsonstruct_geometry.Geometry.nL     = 5;
jsonstruct_geometry.Geometry.nas    = 20;
fjv = flattenJsonStruct(jsonstruct_geometry);

jsonfilename = fullfile('Examples', 'jsondatafiles', 'cc_discharge_control.json');
jsonstruct_control = parseBattmoJson(jsonfilename);

jsonfilename = fullfile('Examples', 'jsondatafiles', 'simulation_parameters.json');
jsonstruct_simparams = parseBattmoJson(jsonfilename);

jsonstruct = mergeJsonStructs({jsonstruct_geometry , ...
    jsonstruct_material , ...
    jsonstruct_control  , ...
    jsonstruct_simparams}, 'warn', false);

model = setupModelFromJson(jsonstruct);

% Define a fun color scheme (custom color scheme)
funColors = [
    0.9, 0.1, 0.5;  % Bright Pink
    0.9, 0.5, 0.1   % Bright Orange
    0.9, 0.5, 0.1   % Bright Orange
    0.2, 0.5, 0.2;  % Forest green
    0.1, 0.3, 0.8;  % Dark Blue
    
    
];

% Ensure the same figure is used for all frames
for i = 1:num_frames
    % Set current figure
    set(0, 'CurrentFigure', figure_handle);  % Ensure we use the same figure
    cla;  % Clear the axes, not the figure
    
    % Use the wrapper function to plot the battery model with your custom colors
    plotBatteryGridWithColors(model, 'fig', figure_handle, 'setstyle', false, 'colors', funColors);
    
    axis tight;
    
    % Remove the legend if it gets created
    legend off;

    % Remove the axes and background
    axis off;  % Hide the axes
    set(gca, 'color', 'none');  % Set the axes background to transparent
    set(gcf, 'color', 'none');  % Set the figure background to transparent
    
    % Rotate the view by 10 degrees per step
    view(45 + 10 * i, 45);
    drawnow;  % Render the updated plot
    
    % Capture the frame
    frame = getframe(figure_handle);  % Capture frame from the same figure
    im = frame2im(frame);  % Convert the frame to an image
    [imind, cm] = rgb2ind(im, 256);  % Convert the frame to indexed image

    % Write the frame to the GIF
    if i == 1
        % First frame: create the GIF with infinite loop
        imwrite(imind, cm, gif_filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.25);
    else
        % Subsequent frames: append to the GIF
        imwrite(imind, cm, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
    end
end




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

    if opt.axisLabels
        xlabel 'x';
        if G.griddim > 1
            ylabel 'y';
            if G.griddim > 2
                zlabel 'z';
            end
        end
    end

    if opt.shortLegendText
        legtext = cellfun(@(s) strrep(s, 'positive electrode', 'PE'), legtext, 'un', false);
        legtext = cellfun(@(s) strrep(s, 'negative electrode', 'NE'), legtext, 'un', false);
    end

    legend off;  % Remove legend to avoid duplication
    axis tight;

    if opt.setstyle
        setFigureStyle('quantity', 'single');
    end

    drawnow();
end

