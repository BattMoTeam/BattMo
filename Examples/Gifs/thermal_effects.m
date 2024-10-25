%{
    This is a script to create a gif showing thermal effects.
%}

% Define the output GIF file name
gif_filename = 'Examples/Gifs/thermal_effects.gif';

jsonstruct_material = parseBattmoJson('Examples/jsondatafiles/sample_input.json');
jsonstruct_material.include_current_collectors = true;

jsonstruct_geometry = parseBattmoJson('Examples/JsonDataFiles/geometry3d.json');
disp(jsonstruct_geometry)

jsonstruct = mergeJsonStructs({jsonstruct_geometry , ...
    jsonstruct_material});

jsonstruct.use_thermal = true;

model = setupModelFromJson(jsonstruct);
output = runBatteryJson(jsonstruct);

%%

% Create the figure once and pass it into plotBatteryGrid
figure_handle = figure('Color', 'none', 'Renderer', 'opengl');  % Create and reuse this figure with no background
% Number of states in the output
num_states = length(output.states)/10;  % Use output.states to get the correct number of states
interval = 1;  % Interval for states to skip

% Loop over every 50th state
for i = 1:interval:num_states
    % Set current figure
    set(0, 'CurrentFigure', figure_handle);  % Ensure we use the same figure
    cla;  % Clear the axes, not the figure
    
    % Access the current state
    state = output.states{i};  % Get the current state

    % Plot temperature data
    plotCellData(model.ThermalModel.grid, ...
        state.ThermalModel.T + T0);

    % Remove the legend if it gets created
    legend off;

    % Remove the axes and background
    axis off;  % Hide the axes
    set(gca, 'color', 'none');  % Set the axes background to transparent
    set(gcf, 'color', 'none');  % Set the figure background to transparent

    title('Temperature / °C');
    axis tight
    view([50, 50]);

    drawnow;  % Update the figure
    
    % Capture the frame and append to the GIF
    frame = getframe(figure_handle);  % Capture the current figure
    im = frame2im(frame);  % Convert the frame to an image
    im = im(:, :, 1:3);  % Remove alpha channel if exists
    [imind, cm] = rgb2ind(im, 256);  % Convert the frame to indexed image
    
    % Write the frame to the GIF
    if i == 1
        % First frame: create the GIF with infinite loop
        imwrite(imind, cm, gif_filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.2);
    else
        % Subsequent frames: append to the GIF
        imwrite(imind, cm, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.2);
    end
end
