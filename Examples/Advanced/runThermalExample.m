%% Example with thermal effects 

%% Setup material properties

jsonfilename = fullfile('ParameterData'        , ...
                        'BatteryCellParameters', ...
                        'LithiumIonBatteryCell', ...
                        'lithium_ion_battery_nmc_graphite.json');
jsonstruct_material = parseBattmoJson(jsonfilename);

jsonstruct_material.include_current_collectors = true;

%% Setup geometry
jsonfilename = fullfile('Examples'     , ...
                        'JsonDataFiles', ...
                        'geometry3d.json');
jsonstruct_geometry = parseBattmoJson(jsonfilename);

%% Setup Control
jsonfilename = fullfile('Examples', 'JsonDataFiles', 'cc_discharge_control.json');
jsonstruct_control = parseBattmoJson(jsonfilename);

jsonstruct = mergeJsonStructs({jsonstruct_geometry , ...
                               jsonstruct_material , ...
                               jsonstruct_control});

%% Plot the extrnal heat transfer coefficient

model = setupModelFromJson(jsonstruct);

% We index of the faces that are coupled to thermally to the exterior
extfaceind = model.ThermalModel.couplingTerm.couplingfaces;
G          = model.ThermalModel.G;
nf         = model.ThermalModel.G.faces.num;

% We create a vector with one value per face with value equal to the external heat transfer coefficient for the external
% face
val = NaN(nf, 1);
val(extfaceind) = model.ThermalModel.externalHeatTransferCoefficient;

figure
plotFaceData(G, val, 'edgecolor', 'black');
axis equal
view([50, 20]);
title('External Heat Transfer Coefficient / W/s/m^2')
colorbar

%% Run the simulation

output = runBatteryJson(jsonstruct);

%% plot voltage

time = output.time;
E    = output.E;

set(0, 'defaulttextfontsize', 15);
set(0, 'defaultaxesfontsize', 15);
set(0, 'defaultlinelinewidth', 3);

figure
plot(time/hour, E)
title('Voltage / V');
xlabel('time / h');
ylabel('voltage / V');


%% Plot the minimum and maximum values of the temperature

Tabs = PhysicalConstants.Tabs;

states = output.states;

Tmin = cellfun(@(state) min(state.ThermalModel.T + Tabs), states);
Tmax = cellfun(@(state) max(state.ThermalModel.T + Tabs), states);

figure
hold on
plot(time / hour, Tmin, 'displayname', 'min T');
plot(time / hour, Tmax, 'displayname', 'max T');
title('Temperature / C')
xlabel('time / h');
ylabel('Temperature / C');

legend show

%% change setup for cooling coefficient

jsonstruct.ThermalModel.externalHeatTransferCoefficientTab = 100;
jsonstruct.ThermalModel.externalHeatTransferCoefficient = 0;

output = runBatteryJson(jsonstruct);

%% Plot the minimum and maximum values of the temperature

states = output.states;

Tmin = cellfun(@(state) min(state.ThermalModel.T + Tabs), states);
Tmax = cellfun(@(state) max(state.ThermalModel.T + Tabs), states);
time = output.time;

figure
hold on
plot(time / hour, Tmin, 'displayname', 'min T');
plot(time / hour, Tmax, 'displayname', 'max T');
title('Temperature / C')
xlabel('time / h');
ylabel('Temperature / C');

legend show

%% plot final temperature disctribution

state = states{end}
figure
plotCellData(model.ThermalModel.G, ...
             state.ThermalModel.T + Tabs);
colorbar
title('Temperature / C');
view([50, 50]);

