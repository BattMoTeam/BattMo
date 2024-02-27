% function [fig] = plotDashboard(model, states, varargin)
% %UNTITLED8 Summary of this function goes here
% %   Detailed explanation goes here

time = output.time/hour;

x = output.model.grid.cells.centroids/(micro*meter);
x_ne = output.model.NegativeElectrode.grid.cells.centroids/(micro*meter);
x_sep = output.model.Separator.grid.cells.centroids/(micro*meter);
x_pe = output.model.PositiveElectrode.grid.cells.centroids/(micro*meter);

% get the centroid locations and the values for the electrolyte
% concentration at the given timestep
% Initialize an empty matrix to store 'c' values
numStates = numel(output.states);  % Get the number of elements in the 'states' array
maxSize = max(cellfun(@(x) numel(x.Electrolyte.c), output.states)); % Get the maximum size of 'c' arrays
c_elyte = NaN(numStates, maxSize);  % Initialize the matrix with NaNs
c_ne = NaN(numStates, maxSize);  % Initialize the matrix with NaNs
c_pe = NaN(numStates, maxSize);  % Initialize the matrix with NaNs
phi_elyte = NaN(numStates, maxSize);  % Initialize the matrix with NaNs
phi_ne = NaN(numStates, maxSize);  % Initialize the matrix with NaNs
phi_pe = NaN(numStates, maxSize);  % Initialize the matrix with NaNs

% Loop through each element in 'states' and extract 'c' values
for i = 1:numStates
    c_elyte(i, 1:numel(output.states{i}.Electrolyte.c)) = output.states{i}.Electrolyte.c;
    c_ne(i, 1:numel(output.states{i}.NegativeElectrode.Coating.ActiveMaterial.SolidDiffusion.cSurface)) = output.states{i}.NegativeElectrode.Coating.ActiveMaterial.SolidDiffusion.cSurface;
    c_pe(i, 1:numel(output.states{i}.PositiveElectrode.Coating.ActiveMaterial.SolidDiffusion.cSurface)) = output.states{i}.PositiveElectrode.Coating.ActiveMaterial.SolidDiffusion.cSurface;
    phi_elyte(i, 1:numel(output.states{i}.Electrolyte.phi)) = output.states{i}.Electrolyte.phi;
    phi_ne(i, 1:numel(output.states{i}.NegativeElectrode.Coating.phi)) = output.states{i}.NegativeElectrode.Coating.phi;
    phi_pe(i, 1:numel(output.states{i}.PositiveElectrode.Coating.phi)) = output.states{i}.PositiveElectrode.Coating.phi;
end

% plot the concentration values at the given grid centroid locations
figure('Units', 'centimeters', 'Position', [0, 0, 29.7, 21]);
contourf(x, time, c_ne ./ jsonstruct.NegativeElectrode.Coating.ActiveMaterial.Interface.saturationConcentration, 20, 'LineStyle', 'none')
xlim([min(x_ne), max(x_ne)])
xlabel('Position  /  µm')
ylabel('Time  /  h')
title('Negative Electrode Lithiation  /  -')
colorbar()
cm = flipud(crameri('lajolla'));
colormap(cm)
set(gca, 'FontSize', 18);

figure('Units', 'centimeters', 'Position', [0, 0, 29.7, 21]);
contourf(x, time, c_elyte, 20, 'LineStyle', 'none')
xlabel('Position  /  µm')
ylabel('Time  /  h')
title('Electrolyte Concentration  /  mol \cdot m^{-3}')
colorbar()
cm = cmocean('curl', 'pivot', 1000);
colormap(cm)
set(gca, 'FontSize', 18);

figure('Units', 'centimeters', 'Position', [0, 0, 29.7, 21]);
contourf(x+max(x_sep), time, c_pe ./ jsonstruct.PositiveElectrode.Coating.ActiveMaterial.Interface.saturationConcentration, 20, 'LineStyle', 'none')
xlim([min(x_pe), max(x_pe)])
xlabel('Position  /  µm')
ylabel('Time  /  h')
title('Positive Electrode Lithiation  /  -')
colorbar()
cm = crameri('nuuk');
colormap(cm)
set(gca, 'FontSize', 18);

% subplot(2,3,4), contourf(x, time, phi_ne, 20, 'LineWidth',0.1)
% xlim([min(x_ne), max(x_ne)])
% xlabel('Position  /  µm')
% ylabel('Time  /  h')
% title('Negative Electrode Potential')
% colorbar()
% 
% subplot(2,3,5), contourf(x, time, phi_elyte, 20, 'LineWidth',0.1)
% xlabel('Position  /  µm')
% ylabel('Time  /  h')
% title('Electrolyte Potential')
% colorbar()
% 
% subplot(2,3,6), contourf(x+max(x_sep), time, phi_pe, 20, 'LineWidth',0.1)
% xlim([min(x_pe), max(x_pe)])
% xlabel('Position  /  µm')
% ylabel('Time  /  h')
% title('Positive Electrode Potential')
% colorbar()


%
%end