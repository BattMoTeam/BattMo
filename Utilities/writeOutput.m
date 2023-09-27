function [] = writeOutput(model, states, name)
%WRITEOUTPUT Summary of this function goes here
%   Detailed explanation goes here

ne      = 'NegativeElectrode';
pe      = 'PositiveElectrode';
elyte   = 'Electrolyte';
thermal = 'ThermalModel';
am      = 'ActiveMaterial';
itf     = 'Interface';
sd      = 'SolidDiffusion';
ctrl    = 'Control';
cc      = 'CurrentCollector';

directory = fileparts(which('startupBattMo'));
file_name = string([directory, '\', name]);

% Check if the file exists, and if so, delete it
if exist(file_name, 'file') == 2
    % Delete the file
    delete(file_name);
end

[~] = H5F.create(file_name);

CRate = model.Control.CRate;

% Instantiate state variables
time = zeros(size(states));
voltage = zeros(size(states));
current = zeros(size(states));

ne_electricPotential = zeros(length(states), length(states{1}.(ne).(am).phi));
ne_concentration = zeros(length(states), length(states{1}.(ne).(am).(sd).c));
ne_surfaceConcentration = zeros(length(states), length(states{1}.(ne).(am).(sd).cSurface));

pe_electricPotential = zeros(length(states), length(states{1}.(pe).(am).phi));
pe_concentration = zeros(length(states), length(states{1}.(pe).(am).(sd).c));
pe_surfaceConcentration = zeros(length(states), length(states{1}.(pe).(am).(sd).cSurface));

elyte_electricPotential = zeros(length(states), length(states{1}.(elyte).phi));
elyte_concentration = zeros(length(states), length(states{1}.(elyte).c));

% Process the states array
for ind = 1:length(states)
    % Get global state variables
    time(ind, 1) = states{ind}.time;
    voltage(ind, 1) = states{ind}.Control.E;
    current(ind, 1) = states{ind}.Control.I;

    % Get local state variables - Negative Electrode
    ne_electricPotential(ind, :) = states{ind}.(ne).(am).phi;
    ne_surfaceConcentration(ind, :) = states{ind}.(ne).(am).(sd).cSurface;
    ne_concentration(ind, :) = states{ind}.(ne).(am).(sd).c;

    % Get local state variables - Positive Electrode
    pe_electricPotential(ind, :) = states{ind}.(pe).(am).phi;
    pe_surfaceConcentration(ind, :) = states{ind}.(pe).(am).(sd).cSurface;
    pe_concentration(ind, :) = states{ind}.(pe).(am).(sd).c;

    % Get local state variables - Electrolyte
    elyte_electricPotential(ind, :) = states{ind}.(elyte).phi;
    elyte_concentration(ind, :) = states{ind}.(elyte).c;

end

% Instantiate global state variables
h5create(file_name, '/global/quantities/time', size(time));
h5create(file_name, '/global/quantities/voltage', size(voltage));
h5create(file_name, '/global/quantities/current', size(current));

% Instantiate local state variables
h5create(file_name, '/local/NegativeElectrode/quantities/surfaceConcentration', size(ne_surfaceConcentration));
h5create(file_name, '/local/NegativeElectrode/quantities/concentration', size(ne_concentration));
h5create(file_name, '/local/NegativeElectrode/quantities/electricPotential', size(ne_electricPotential));
h5create(file_name, '/local/PositiveElectrode/quantities/surfaceConcentration', size(pe_surfaceConcentration));
h5create(file_name, '/local/PositiveElectrode/quantities/concentration', size(pe_concentration));
h5create(file_name, '/local/PositiveElectrode/quantities/electricPotential', size(pe_electricPotential));
h5create(file_name, '/local/Electrolyte/quantities/concentration', size(elyte_concentration));
h5create(file_name, '/local/Electrolyte/quantities/electricPotential', size(elyte_electricPotential));

% Instantiate specifications
h5create(file_name, '/global/specification/cRate', size(CRate));

% Instantiate global grid properties
h5create(file_name, '/global/grid/cells/num', size(model.G.cells.num));
h5create(file_name, '/global/grid/cells/facePos', size(model.G.cells.facePos));
h5create(file_name, '/global/grid/cells/indexMap', size(model.G.cells.indexMap));
h5create(file_name, '/global/grid/cells/faces', size(model.G.cells.faces));
h5create(file_name, '/global/grid/cells/volumes', size(model.G.cells.volumes));
h5create(file_name, '/global/grid/cells/centroids', size(model.G.cells.centroids));
h5create(file_name, '/global/grid/faces/num', size(model.G.faces.num));
h5create(file_name, '/global/grid/faces/nodePos', size(model.G.faces.nodePos));
h5create(file_name, '/global/grid/faces/neighbors', size(model.G.faces.neighbors));
h5create(file_name, '/global/grid/faces/nodes', size(model.G.faces.nodes));
h5create(file_name, '/global/grid/faces/areas', size(model.G.faces.areas));
h5create(file_name, '/global/grid/faces/normals', size(model.G.faces.normals));
h5create(file_name, '/global/grid/faces/centroids', size(model.G.faces.centroids));
h5create(file_name, '/global/grid/nodes/num', size(model.G.nodes.num));
h5create(file_name, '/global/grid/nodes/coords', size(model.G.nodes.coords));
h5create(file_name, '/global/grid/cartDims', size(model.G.cartDims));
h5create(file_name, '/global/grid/griddim', size(model.G.griddim));

% Instantiate negative electrode grid properties
h5create(file_name, '/local/NegativeElectrode/grid/cells/num', size(model.(ne).G.cells.num));
h5create(file_name, '/local/NegativeElectrode/grid/cells/facePos', size(model.(ne).G.cells.facePos));
h5create(file_name, '/local/NegativeElectrode/grid/cells/indexMap', size(model.(ne).G.cells.indexMap));
h5create(file_name, '/local/NegativeElectrode/grid/cells/faces', size(model.(ne).G.cells.faces));
h5create(file_name, '/local/NegativeElectrode/grid/cells/volumes', size(model.(ne).G.cells.volumes));
h5create(file_name, '/local/NegativeElectrode/grid/cells/centroids', size(model.(ne).G.cells.centroids));
h5create(file_name, '/local/NegativeElectrode/grid/faces/num', size(model.(ne).G.faces.num));
h5create(file_name, '/local/NegativeElectrode/grid/faces/nodePos', size(model.(ne).G.faces.nodePos));
h5create(file_name, '/local/NegativeElectrode/grid/faces/neighbors', size(model.(ne).G.faces.neighbors));
h5create(file_name, '/local/NegativeElectrode/grid/faces/nodes', size(model.(ne).G.faces.nodes));
h5create(file_name, '/local/NegativeElectrode/grid/faces/areas', size(model.(ne).G.faces.areas));
h5create(file_name, '/local/NegativeElectrode/grid/faces/normals', size(model.(ne).G.faces.normals));
h5create(file_name, '/local/NegativeElectrode/grid/faces/centroids', size(model.(ne).G.faces.centroids));
h5create(file_name, '/local/NegativeElectrode/grid/nodes/num', size(model.(ne).G.nodes.num));
h5create(file_name, '/local/NegativeElectrode/grid/nodes/coords', size(model.(ne).G.nodes.coords));
h5create(file_name, '/local/NegativeElectrode/grid/cartDims', size(model.(ne).G.cartDims));
h5create(file_name, '/local/NegativeElectrode/grid/griddim', size(model.(ne).G.griddim));

% Instantiate positive electrode grid properties
h5create(file_name, '/local/PositiveElectrode/grid/cells/num', size(model.(pe).G.cells.num));
h5create(file_name, '/local/PositiveElectrode/grid/cells/facePos', size(model.(pe).G.cells.facePos));
h5create(file_name, '/local/PositiveElectrode/grid/cells/indexMap', size(model.(pe).G.cells.indexMap));
h5create(file_name, '/local/PositiveElectrode/grid/cells/faces', size(model.(pe).G.cells.faces));
h5create(file_name, '/local/PositiveElectrode/grid/cells/volumes', size(model.(pe).G.cells.volumes));
h5create(file_name, '/local/PositiveElectrode/grid/cells/centroids', size(model.(pe).G.cells.centroids));
h5create(file_name, '/local/PositiveElectrode/grid/faces/num', size(model.(pe).G.faces.num));
h5create(file_name, '/local/PositiveElectrode/grid/faces/nodePos', size(model.(pe).G.faces.nodePos));
h5create(file_name, '/local/PositiveElectrode/grid/faces/neighbors', size(model.(pe).G.faces.neighbors));
h5create(file_name, '/local/PositiveElectrode/grid/faces/nodes', size(model.(pe).G.faces.nodes));
h5create(file_name, '/local/PositiveElectrode/grid/faces/areas', size(model.(pe).G.faces.areas));
h5create(file_name, '/local/PositiveElectrode/grid/faces/normals', size(model.(pe).G.faces.normals));
h5create(file_name, '/local/PositiveElectrode/grid/faces/centroids', size(model.(pe).G.faces.centroids));
h5create(file_name, '/local/PositiveElectrode/grid/nodes/num', size(model.(pe).G.nodes.num));
h5create(file_name, '/local/PositiveElectrode/grid/nodes/coords', size(model.(pe).G.nodes.coords));
h5create(file_name, '/local/PositiveElectrode/grid/cartDims', size(model.(pe).G.cartDims));
h5create(file_name, '/local/PositiveElectrode/grid/griddim', size(model.(pe).G.griddim));

% Instantiate electrolyte grid properties
h5create(file_name, '/local/Electrolyte/grid/cells/num', size(model.(elyte).G.cells.num));
h5create(file_name, '/local/Electrolyte/grid/cells/facePos', size(model.(elyte).G.cells.facePos));
h5create(file_name, '/local/Electrolyte/grid/cells/indexMap', size(model.(elyte).G.cells.indexMap));
h5create(file_name, '/local/Electrolyte/grid/cells/faces', size(model.(elyte).G.cells.faces));
h5create(file_name, '/local/Electrolyte/grid/cells/volumes', size(model.(elyte).G.cells.volumes));
h5create(file_name, '/local/Electrolyte/grid/cells/centroids', size(model.(elyte).G.cells.centroids));
h5create(file_name, '/local/Electrolyte/grid/faces/num', size(model.(elyte).G.faces.num));
h5create(file_name, '/local/Electrolyte/grid/faces/nodePos', size(model.(elyte).G.faces.nodePos));
h5create(file_name, '/local/Electrolyte/grid/faces/neighbors', size(model.(elyte).G.faces.neighbors));
h5create(file_name, '/local/Electrolyte/grid/faces/nodes', size(model.(elyte).G.faces.nodes));
h5create(file_name, '/local/Electrolyte/grid/faces/areas', size(model.(elyte).G.faces.areas));
h5create(file_name, '/local/Electrolyte/grid/faces/normals', size(model.(elyte).G.faces.normals));
h5create(file_name, '/local/Electrolyte/grid/faces/centroids', size(model.(elyte).G.faces.centroids));
h5create(file_name, '/local/Electrolyte/grid/nodes/num', size(model.(elyte).G.nodes.num));
h5create(file_name, '/local/Electrolyte/grid/nodes/coords', size(model.(elyte).G.nodes.coords));
h5create(file_name, '/local/Electrolyte/grid/cartDims', size(model.(elyte).G.cartDims));
h5create(file_name, '/local/Electrolyte/grid/griddim', size(model.(elyte).G.griddim));



% Write global state variables
h5write(file_name, '/global/quantities/time', time);
h5write(file_name, '/global/quantities/voltage', voltage);
h5write(file_name, '/global/quantities/current', current);

% Write local state variables
h5write(file_name, '/local/NegativeElectrode/quantities/surfaceConcentration', ne_surfaceConcentration);
h5write(file_name, '/local/NegativeElectrode/quantities/concentration', ne_concentration);
h5write(file_name, '/local/NegativeElectrode/quantities/electricPotential', ne_electricPotential);
h5write(file_name, '/local/PositiveElectrode/quantities/surfaceConcentration', pe_surfaceConcentration);
h5write(file_name, '/local/PositiveElectrode/quantities/concentration', pe_concentration);
h5write(file_name, '/local/PositiveElectrode/quantities/electricPotential', pe_electricPotential);
h5write(file_name, '/local/Electrolyte/quantities/concentration', elyte_concentration);
h5write(file_name, '/local/Electrolyte/quantities/electricPotential', elyte_electricPotential);

% Write specifications
h5write(file_name, '/global/specification/cRate', CRate);

% Write global grid properties
h5write(file_name, '/global/grid/cells/num', model.G.cells.num);
h5write(file_name, '/global/grid/cells/facePos', model.G.cells.facePos);
h5write(file_name, '/global/grid/cells/indexMap', model.G.cells.indexMap);
h5write(file_name, '/global/grid/cells/faces', model.G.cells.faces);
h5write(file_name, '/global/grid/cells/volumes', model.G.cells.volumes);
h5write(file_name, '/global/grid/cells/centroids', model.G.cells.centroids);
h5write(file_name, '/global/grid/faces/num', model.G.faces.num);
h5write(file_name, '/global/grid/faces/nodePos', model.G.faces.nodePos);
h5write(file_name, '/global/grid/faces/neighbors', model.G.faces.neighbors);
h5write(file_name, '/global/grid/faces/nodes', model.G.faces.nodes);
h5write(file_name, '/global/grid/faces/areas', model.G.faces.areas);
h5write(file_name, '/global/grid/faces/normals', model.G.faces.normals);
h5write(file_name, '/global/grid/faces/centroids', model.G.faces.centroids);
h5write(file_name, '/global/grid/nodes/num', model.G.nodes.num);
h5write(file_name, '/global/grid/nodes/coords', model.G.nodes.coords);
h5write(file_name, '/global/grid/cartDims', model.G.cartDims);
h5write(file_name, '/global/grid/griddim', model.G.griddim);

% Write negative electrode grid properties
h5write(file_name, '/local/NegativeElectrode/grid/cells/num', model.(ne).G.cells.num);
h5write(file_name, '/local/NegativeElectrode/grid/cells/facePos', model.(ne).G.cells.facePos);
h5write(file_name, '/local/NegativeElectrode/grid/cells/indexMap', model.(ne).G.cells.indexMap);
h5write(file_name, '/local/NegativeElectrode/grid/cells/faces', model.(ne).G.cells.faces);
h5write(file_name, '/local/NegativeElectrode/grid/cells/volumes', model.(ne).G.cells.volumes);
h5write(file_name, '/local/NegativeElectrode/grid/cells/centroids', model.(ne).G.cells.centroids);
h5write(file_name, '/local/NegativeElectrode/grid/faces/num', model.(ne).G.faces.num);
h5write(file_name, '/local/NegativeElectrode/grid/faces/nodePos', model.(ne).G.faces.nodePos);
h5write(file_name, '/local/NegativeElectrode/grid/faces/neighbors', model.(ne).G.faces.neighbors);
h5write(file_name, '/local/NegativeElectrode/grid/faces/nodes', model.(ne).G.faces.nodes);
h5write(file_name, '/local/NegativeElectrode/grid/faces/areas', model.(ne).G.faces.areas);
h5write(file_name, '/local/NegativeElectrode/grid/faces/normals', model.(ne).G.faces.normals);
h5write(file_name, '/local/NegativeElectrode/grid/faces/centroids', model.(ne).G.faces.centroids);
h5write(file_name, '/local/NegativeElectrode/grid/nodes/num', model.(ne).G.nodes.num);
h5write(file_name, '/local/NegativeElectrode/grid/nodes/coords', model.(ne).G.nodes.coords);
h5write(file_name, '/local/NegativeElectrode/grid/cartDims', model.(ne).G.cartDims);
h5write(file_name, '/local/NegativeElectrode/grid/griddim', model.(ne).G.griddim);

% Write positive electrode grid properties
h5write(file_name, '/local/PositiveElectrode/grid/cells/num', model.(pe).G.cells.num);
h5write(file_name, '/local/PositiveElectrode/grid/cells/facePos', model.(pe).G.cells.facePos);
h5write(file_name, '/local/PositiveElectrode/grid/cells/indexMap', model.(pe).G.cells.indexMap);
h5write(file_name, '/local/PositiveElectrode/grid/cells/faces', model.(pe).G.cells.faces);
h5write(file_name, '/local/PositiveElectrode/grid/cells/volumes', model.(pe).G.cells.volumes);
h5write(file_name, '/local/PositiveElectrode/grid/cells/centroids', model.(pe).G.cells.centroids);
h5write(file_name, '/local/PositiveElectrode/grid/faces/num', model.(pe).G.faces.num);
h5write(file_name, '/local/PositiveElectrode/grid/faces/nodePos', model.(pe).G.faces.nodePos);
h5write(file_name, '/local/PositiveElectrode/grid/faces/neighbors', model.(pe).G.faces.neighbors);
h5write(file_name, '/local/PositiveElectrode/grid/faces/nodes', model.(pe).G.faces.nodes);
h5write(file_name, '/local/PositiveElectrode/grid/faces/areas', model.(pe).G.faces.areas);
h5write(file_name, '/local/PositiveElectrode/grid/faces/normals', model.(pe).G.faces.normals);
h5write(file_name, '/local/PositiveElectrode/grid/faces/centroids', model.(pe).G.faces.centroids);
h5write(file_name, '/local/PositiveElectrode/grid/nodes/num', model.(pe).G.nodes.num);
h5write(file_name, '/local/PositiveElectrode/grid/nodes/coords', model.(pe).G.nodes.coords);
h5write(file_name, '/local/PositiveElectrode/grid/cartDims', model.(pe).G.cartDims);
h5write(file_name, '/local/PositiveElectrode/grid/griddim', model.(pe).G.griddim);

% Write electrolyte grid properties
h5write(file_name, '/local/Electrolyte/grid/cells/num', model.(elyte).G.cells.num);
h5write(file_name, '/local/Electrolyte/grid/cells/facePos', model.(elyte).G.cells.facePos);
h5write(file_name, '/local/Electrolyte/grid/cells/indexMap', model.(elyte).G.cells.indexMap);
h5write(file_name, '/local/Electrolyte/grid/cells/faces', model.(elyte).G.cells.faces);
h5write(file_name, '/local/Electrolyte/grid/cells/volumes', model.(elyte).G.cells.volumes);
h5write(file_name, '/local/Electrolyte/grid/cells/centroids', model.(elyte).G.cells.centroids);
h5write(file_name, '/local/Electrolyte/grid/faces/num', model.(elyte).G.faces.num);
h5write(file_name, '/local/Electrolyte/grid/faces/nodePos', model.(elyte).G.faces.nodePos);
h5write(file_name, '/local/Electrolyte/grid/faces/neighbors', model.(elyte).G.faces.neighbors);
h5write(file_name, '/local/Electrolyte/grid/faces/nodes', model.(elyte).G.faces.nodes);
h5write(file_name, '/local/Electrolyte/grid/faces/areas', model.(elyte).G.faces.areas);
h5write(file_name, '/local/Electrolyte/grid/faces/normals', model.(elyte).G.faces.normals);
h5write(file_name, '/local/Electrolyte/grid/faces/centroids', model.(elyte).G.faces.centroids);
h5write(file_name, '/local/Electrolyte/grid/nodes/num', model.(elyte).G.nodes.num);
h5write(file_name, '/local/Electrolyte/grid/nodes/coords', model.(elyte).G.nodes.coords);
h5write(file_name, '/local/Electrolyte/grid/cartDims', model.(elyte).G.cartDims);
h5write(file_name, '/local/Electrolyte/grid/griddim', model.(elyte).G.griddim);


% if exist(specs_name, 'file') == 2
%     % Delete the file
%     specs = load(specs_name, 'specs');
%     elde_specs = density_computation();
% 
%     cell_mass = specs.specs.mass .* 1000;
%     ne_mass = (specs.specs.masses.NegativeElectrode.ActiveMaterial.val + specs.specs.masses.NegativeElectrode.CurrentCollector.val) .* 1000;
%     pe_mass = (specs.specs.masses.PositiveElectrode.ActiveMaterial.val + specs.specs.masses.PositiveElectrode.CurrentCollector.val) .* 1000;
%     elyte_mass = (specs.specs.masses.Electrolyte.val) .* 1000;
%     sep_mass = (specs.specs.masses.Electrolyte.Separator.val) .* 1000;
%     np_ratio = specs.specs.NPratio;
% 
%     ne_spec_cap = specs.specs.cap_neg/3600/specs.specs.masses.NegativeElectrode.ActiveMaterial.val;
%     pe_spec_cap = specs.specs.cap_pos/3600/specs.specs.masses.PositiveElectrode.ActiveMaterial.val;
% 
%     cell_capacity = specs.specs.cap ./ 3600;
%     cell_energy = specs.specs.effectiveEnergy ./ 3600;
%     cell_specificEnergy = cell_energy ./ (cell_mass ./ 1000);
% 
%     ne_thickness = output_data.output.jsonstruct.NegativeElectrode.ActiveMaterial.thickness;
%     pe_thickness = output_data.output.jsonstruct.PositiveElectrode.ActiveMaterial.thickness;
% 
%     ne_am_density = elde_specs.densities.NegativeElectrode.ActiveMaterial;
%     ne_bind_density = elde_specs.densities.NegativeElectrode.Binder;
%     ne_add_density = elde_specs.densities.NegativeElectrode.ConductingAdditive;
% 
%     ne_am_wt = elde_specs.compositions.NegativeElectrode.ActiveMaterial;
%     ne_bind_wt = elde_specs.compositions.NegativeElectrode.Binder;
%     ne_add_wt = elde_specs.compositions.NegativeElectrode.ConductingAdditive;
% 
%     pe_am_density = elde_specs.densities.PositiveElectrode.ActiveMaterial;
%     pe_bind_density = elde_specs.densities.PositiveElectrode.Binder;
%     pe_add_density = elde_specs.densities.PositiveElectrode.ConductingAdditive;
% 
%     pe_am_wt = elde_specs.compositions.PositiveElectrode.ActiveMaterial;
%     pe_bind_wt = elde_specs.compositions.PositiveElectrode.Binder;
%     pe_add_wt = elde_specs.compositions.PositiveElectrode.ConductingAdditive;
% 
%     ne_coating_density = elde_specs.density.NegativeElectrode ./ 1000;
%     pe_coating_density = elde_specs.density.PositiveElectrode ./ 1000;
% 
%     ne_mass_loading = elde_specs.convertionfactor.NegativeElectrode .* ne_thickness ./ 0.01;
%     pe_mass_loading = elde_specs.convertionfactor.PositiveElectrode .* pe_thickness ./ 0.01;
% 
% 
%     h5create(file_name, '/cell_mass', size(cell_mass));
%     h5create(file_name, '/ne_mass', size(ne_mass));
%     h5create(file_name, '/pe_mass', size(pe_mass));
%     h5create(file_name, '/elyte_mass', size(elyte_mass));
%     h5create(file_name, '/sep_mass', size(sep_mass));
%     h5create(file_name, '/np_ratio', size(np_ratio));
%     h5create(file_name, '/ne_spec_cap', size(ne_spec_cap));
%     h5create(file_name, '/pe_spec_cap', size(pe_spec_cap));
%     h5create(file_name, '/cell_capacity', size(cell_capacity));
%     h5create(file_name, '/cell_energy', size(cell_energy));
%     h5create(file_name, '/cell_specificEnergy', size(cell_specificEnergy));
%     h5create(file_name, '/ne_thickness', size(ne_thickness));
%     h5create(file_name, '/pe_thickness', size(pe_thickness));
%     h5create(file_name, '/ne_am_density', size(ne_am_density));
%     h5create(file_name, '/ne_bind_density', size(ne_bind_density));
%     h5create(file_name, '/ne_add_density', size(ne_add_density));
%     h5create(file_name, '/ne_am_wt', size(ne_am_wt));
%     h5create(file_name, '/ne_bind_wt', size(ne_bind_wt));
%     h5create(file_name, '/ne_add_wt', size(ne_add_wt));
%     h5create(file_name, '/pe_am_density', size(pe_am_density));
%     h5create(file_name, '/pe_bind_density', size(pe_bind_density));
%     h5create(file_name, '/pe_add_density', size(pe_add_density));
%     h5create(file_name, '/pe_am_wt', size(pe_am_wt));
%     h5create(file_name, '/pe_bind_wt', size(pe_bind_wt));
%     h5create(file_name, '/pe_add_wt', size(pe_add_wt));
%     h5create(file_name, '/ne_coating_density', size(ne_coating_density));
%     h5create(file_name, '/pe_coating_density', size(pe_coating_density));
%     h5create(file_name, '/ne_mass_loading', size(ne_mass_loading));
%     h5create(file_name, '/pe_mass_loading', size(pe_mass_loading));
% 
%     h5write(file_name, '/cell_mass', (cell_mass));
%     h5write(file_name, '/ne_mass', (ne_mass));
%     h5write(file_name, '/pe_mass', (pe_mass));
%     h5write(file_name, '/elyte_mass', (elyte_mass));
%     h5write(file_name, '/sep_mass', (sep_mass));
%     h5write(file_name, '/np_ratio', (np_ratio));
%     h5write(file_name, '/ne_spec_cap', (ne_spec_cap));
%     h5write(file_name, '/pe_spec_cap', (pe_spec_cap));
%     h5write(file_name, '/cell_capacity', (cell_capacity));
%     h5write(file_name, '/cell_energy', (cell_energy));
%     h5write(file_name, '/cell_specificEnergy', (cell_specificEnergy));
%     h5write(file_name, '/ne_thickness', (ne_thickness));
%     h5write(file_name, '/pe_thickness', (pe_thickness));
%     h5write(file_name, '/ne_am_density', (ne_am_density));
%     h5write(file_name, '/ne_bind_density', (ne_bind_density));
%     h5write(file_name, '/ne_add_density', (ne_add_density));
%     h5write(file_name, '/ne_am_wt', (ne_am_wt));
%     h5write(file_name, '/ne_bind_wt', (ne_bind_wt));
%     h5write(file_name, '/ne_add_wt', (ne_add_wt));
%     h5write(file_name, '/pe_am_density', (pe_am_density));
%     h5write(file_name, '/pe_bind_density', (pe_bind_density));
%     h5write(file_name, '/pe_add_density', (pe_add_density));
%     h5write(file_name, '/pe_am_wt', (pe_am_wt));
%     h5write(file_name, '/pe_bind_wt', (pe_bind_wt));
%     h5write(file_name, '/pe_add_wt', (pe_add_wt));
%     h5write(file_name, '/ne_coating_density', (ne_coating_density));
%     h5write(file_name, '/pe_coating_density', (pe_coating_density));
%     h5write(file_name, '/ne_mass_loading', (ne_mass_loading));
%     h5write(file_name, '/pe_mass_loading', (pe_mass_loading));


end

