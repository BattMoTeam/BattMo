%% Comparison between coupled and uncoupled thermal simulation
% To compute the thermal response of a battery cell, one can either run a fully 
% coupled simulation or run first an isothermal simulation and then use the output 
% to compute the thermal response in a thermal-only simulation.
% 
% In the second case, we do not get any feedback from the thermal model to the 
% electrochemical model.
% 
% In this example, we compare the solutions obtained both for the voltage and 
% for the temperature on a P4D model.
%% setup material property input
% We use a lithium-ion battery cell with NMC cathode and graphite anode

jsonfilename = fullfile('ParameterData'        , ...
                        'BatteryCellParameters', ...
                        'LithiumIonBatteryCell', ...
                        'lithium_ion_battery_nmc_graphite.json');
jsonstruct_material = parseBattmoJson(jsonfilename);
%% Setup geometry input
% We use a simple 3d-geometry (see image below) with only one layer

jsonfilename = fullfile('Examples'     , ...
                        'JsonDataFiles', ...
                        'geometry3d.json');
jsonstruct_geometry = parseBattmoJson(jsonfilename);
%% Setup Control input
% We use a single discharge scenario

jsonfilename = fullfile('Examples', 'JsonDataFiles', 'cc_discharge_control.json');
jsonstruct_control = parseBattmoJson(jsonfilename);
%% 
% We merge the structures into a single input structure

jsonstruct = mergeJsonStructs({jsonstruct_geometry , ...
                               jsonstruct_material , ...
                               jsonstruct_control}, 'warn', false);
%% Modify input structure
% We modify the input parameters manually to obtain higher temperature increases 
% so that the difference between the two simulation approaches get enhanced.
% 
% We change the heat transfer coefficient to low values

jsonstruct.ThermalModel.externalHeatTransferCoefficientTab = 1e-1;
jsonstruct.ThermalModel.externalHeatTransferCoefficient    = 1e-1;
%% 
% We reduce by a constant factor the specific heat capacities of all the materials

coef = 5e-2;
locations = {{'NegativeElectrode', 'CurrentCollector', 'specificHeatCapacity'}, ...
             {'NegativeElectrode', 'Coating', 'ActiveMaterial', 'specificHeatCapacity'}, ...
             {'PositiveElectrode', 'CurrentCollector', 'specificHeatCapacity'}, ...
             {'PositiveElectrode', 'Coating', 'ActiveMaterial', 'specificHeatCapacity'}, ...
             {'Electrolyte', 'specificHeatCapacity'}};

for iloc = 1 : numel(locations)
    loc = locations{iloc};
    val = getJsonStructField(jsonstruct, loc);
    jsonstruct = setJsonStructField(jsonstruct, loc, coef*val, 'handleMisMatch', 'quiet');
end
%% 
% We change the lower cutoff voltage

jsonstruct.Control.lowerCutoffVoltage = 3.6;
%% Run thermal simulation

output_fullycoupled = runBatteryJson(jsonstruct);
%% Plot mesh
% Now the model is setup, we can plot the grid mesh

plotBatteryGrid(output_fullycoupled.model, 'shortLegendText', true, 'figure', 1);
%% Run iso-thermal simulation
% We modify the input structure and switch to iso-thermal

jsonstruct_isothermal = setJsonStructField(jsonstruct, 'use_thermal', false, 'handleMisMatch', 'quiet');
%% 
% We run the iso-thermal simulation

output_isothermal = runBatteryJson(jsonstruct_isothermal);
%% Setup thermal-only simulation
% We compute the source terms from the output states obtained from the isothermal 
% simulation. The source terms depends on the norm of the charge and mass fluxes 
% and on the reaction rates.

states = output_isothermal.states;
for istate = 1 : numel(states)
    % The functions to compute the source terms are available in the fully coupled model.
    states{istate} = output_fullycoupled.model.evalVarName(states{istate}, {'ThermalModel', 'jHeatSource'});
    states{istate} = output_fullycoupled.model.evalVarName(states{istate}, {'ThermalModel', 'jHeatOhmSource'});
    states{istate} = output_fullycoupled.model.evalVarName(states{istate}, {'ThermalModel', 'jHeatChemicalSource'});
    states{istate} = output_fullycoupled.model.updateThermalIrreversibleReactionSourceTerms(states{istate});
    states{istate} = output_fullycoupled.model.updateThermalReversibleReactionSourceTerms(states{istate});
end
states_heat = states;
%% Setup heat source term
% We use the values of thermal source that we just computed to setup an helper 
% structure which is going to be used to send the heat source to the thermal-only 
% simulation

sourceTerms = cellfun(@(state) state.ThermalModel.jHeatSource, states_heat, 'uniformoutput', false);
disp(size(sourceTerms))

disp(size(states_heat{1}.ThermalModel.jHeatSource))


hss = HeatSourceSetup(sourceTerms, output_isothermal.time);
%% Setup thermal-only model
% The thermal-only model uses the same code as the thermal component model used 
% in the fully coupled simulation

inputparams_thermal = output_fullycoupled.inputparams.ThermalModel;
%% 
% The effective thermal properties depends on the intrinsic thermal property 
% of each of the constituant of the battery. The effective thermal properties 
% are setup when the fully-coupled model is instantiated. We recover those property 
% to use them for the thermal-only model

inputparams_thermal.effectiveThermalConductivity    = output_fullycoupled.model.ThermalModel.effectiveThermalConductivity;
inputparams_thermal.effectiveVolumetricHeatCapacity = output_fullycoupled.model.ThermalModel.effectiveVolumetricHeatCapacity;
effectiveThermalConductivity = inputparams_thermal.effectiveThermalConductivity;
effectiveVolumetricHeatCapacity = inputparams_thermal.effectiveVolumetricHeatCapacity;



disp(inputparams_thermal)
model_thermal = ThermalComponent(inputparams_thermal);

%% Boundary arrays used by thermal-only model (for Julia parity checks).
thermalBoundaryNeighbors = [];
thermalBoundaryAreas = [];
try
    coup = model_thermal.couplingTerm;
    coupcells = coup.couplingcells;
    coupfaces = coup.couplingfaces;

    if ~isempty(coupcells)
        thermalBoundaryNeighbors = coupcells(:);

        if ~isempty(coupfaces)
            allFaceAreas = model_thermal.G.getFaceAreas();
            thermalBoundaryAreas = allFaceAreas(coupfaces(:));
        else
            % 1D volumetric-cooling fallback
            v = model_thermal.G.getVolumes();
            thermalBoundaryAreas = v(coupcells(:));
        end

        fprintf('Extracted thermal boundary arrays from couplingTerm: %d entries.\n', numel(thermalBoundaryNeighbors));
    else
        warning('Thermal couplingTerm.couplingcells is empty.');
    end
catch me
    warning('Thermal boundary array extraction failed: %s', me.message);
end
%% 
% The thermal model is used as the main simulation model. We equip it for simulation

model_thermal.isRootSimulationModel = true;
model_thermal = model_thermal.equipModelForComputation();

%% Setup the simulation schedule
% The simulation schedule contains also the heat source term, which is an external 
% source seen from the thermal-only model.
% 
% We use the same time steps as for the iso-thermal simulation.

times = hss.times;

clear step
step.val     = [times(1); diff(times)];
step.control = ones(numel(times), 1);

clear control
control.src = @(time) hss.eval(time);

schedule = struct('control', control, ...
                  'step'   , step);
%% setup initial state
% The initial state is a constant temperature. For convenience, we just retrieve 
% it from the fully-coupled simulation

initstate = output_fullycoupled.model.setupInitialState(jsonstruct);

clear state0;
state0.T = initstate.ThermalModel.T;
%% run thermal-only simulation

simInput = struct('model'    , model_thermal, ...
                  'initstate', state0, ...
                  'schedule' , schedule);

simsetup = SimulationSetup(simInput);

states_thermal = simsetup.run();
%% Plotting
% We plot the evolution of the maximum temperature for the full-coupled and 
% thermal-only simulations and the discharge voltage curves for the iso-thermal 
% and fully-coupled simulations.

T0 = PhysicalConstants.absoluteTemperature;

time   = output_fullycoupled.time;
E      = output_fullycoupled.E;
states = output_fullycoupled.states;

Tmax = cellfun(@(state) max(state.ThermalModel.T + T0), states);

figure(2)
hold on
plot(time / hour, Tmax, 'displayname', 'max T (fully coupled)');
title('Temperature / C')
xlabel('time / h');
ylabel('Temperature / C');

figure(3)
hold on
plot(time/hour, E, 'displayname', 'fully coupled')
title('Voltage / V');
xlabel('time / h');
ylabel('voltage / V');

states = states_thermal;

time = cellfun(@(state) state.time, output_isothermal.states);
E    = cellfun(@(state) state.Control.E, output_isothermal.states);
Tmin = cellfun(@(state) min(state.T + T0), states);
Tmax = cellfun(@(state) max(state.T + T0), states);

figure(2)
plot(time / hour, Tmax, 'displayname', 'max T (decoupled)');
title('Temperature / C')
xlabel('time / h');
ylabel('Temperature / C');

legend show

figure(3)
plot(time/hour, E, 'displayname', 'isothermal')
title('Voltage / V');
xlabel('time / h');
ylabel('voltage / V');

legend show

%% Heat-source decomposition plot (same structure as Julia post-processing)
labels = {'Ohmic current collector (pos)', ...
          'Ohmic current collector (neg)', ...
          'Ohmic active material (pos)', ...
          'Ohmic active material (neg)', ...
          'Ohmic (elyte)', ...
          'Reaction reversible (pos)', ...
          'Reaction reversible (neg)', ...
          'Reaction irreversible (pos)', ...
          'Reaction irreversible (neg)', ...
          'Diffusion (elyte)'};

colors = [0.86 0.37 0.34; ...
          0.95 0.62 0.29; ...
          0.99 0.82 0.32; ...
          0.63 0.79 0.36; ...
          0.27 0.73 0.58; ...
          0.30 0.67 0.88; ...
          0.35 0.49 0.89; ...
          0.55 0.41 0.82; ...
          0.86 0.44 0.67; ...
          0.62 0.54 0.47];

nlabels = numel(labels);
nt = numel(states_heat);
P = zeros(nlabels, nt);
Ptot = zeros(1, nt);

model = output_fullycoupled.model;
ne = 'NegativeElectrode'; pe = 'PositiveElectrode';
co = 'Coating'; cc = 'CurrentCollector';
am = 'ActiveMaterial'; sd = 'SolidDiffusion'; itf = 'Interface';
elyte = 'Electrolyte';

for i = 1 : nt
    st = states_heat{i};

    % Total
    Ptot(i) = sum(st.ThermalModel.jHeatSource);

    % Ohmic current collectors
    if model.include_current_collectors
        cc_model = model.(pe).(cc);
        cc_j     = st.(pe).(cc).jFace;
        cc_econd = cc_model.effectiveElectronicConductivity;
        cc_vols  = cc_model.G.getVolumes();
        cc_jsq   = cc_model.G.getCellFluxNorm(cc_j);
        P(1, i)  = sum(cc_vols .* cc_jsq ./ cc_econd);

        cc_model = model.(ne).(cc);
        cc_j     = st.(ne).(cc).jFace;
        cc_econd = cc_model.effectiveElectronicConductivity;
        cc_vols  = cc_model.G.getVolumes();
        cc_jsq   = cc_model.G.getCellFluxNorm(cc_j);
        P(2, i)  = sum(cc_vols .* cc_jsq ./ cc_econd);
    end

    % Ohmic active material (coatings)
    co_model = model.(pe).(co);
    co_j     = st.(pe).(co).jFace;
    co_econd = co_model.effectiveElectronicConductivity;
    co_vols  = co_model.G.getVolumes();
    co_jsq   = co_model.G.getCellFluxNorm(co_j);
    P(3, i)  = sum(co_vols .* co_jsq ./ co_econd);

    co_model = model.(ne).(co);
    co_j     = st.(ne).(co).jFace;
    co_econd = co_model.effectiveElectronicConductivity;
    co_vols  = co_model.G.getVolumes();
    co_jsq   = co_model.G.getCellFluxNorm(co_j);
    P(4, i)  = sum(co_vols .* co_jsq ./ co_econd);

    % Ohmic electrolyte
    el_model  = model.(elyte);
    el_j      = st.(elyte).jFace;
    el_vf     = el_model.volumeFraction;
    el_brug   = el_model.bruggemanCoefficient;
    el_cond   = st.(elyte).conductivity;
    el_econd  = el_cond .* el_vf .^ el_brug;
    el_vols   = el_model.G.getVolumes();
    el_jsq    = el_model.G.getCellFluxNorm(el_j);
    P(5, i)   = sum(el_vols .* el_jsq ./ el_econd);

    % Reaction reversible / irreversible (positive and negative)
    T = st.ThermalModel.T;
    for k = 1 : 2
        if k == 1
            elde = pe;
            irev = 8;
            irevb = 6;
        else
            elde = ne;
            irev = 9;
            irevb = 7;
        end
        F      = model.(elde).(co).(am).(itf).constants.F;
        n      = model.(elde).(co).(am).(itf).numberOfElectronsTransferred;
        cmap   = model.(elde).(co).G.mappings.cellmap;
        vols   = model.(elde).(co).G.getVolumes();
        Rvol   = st.(elde).(co).(am).(sd).Rvol;
        eta    = st.(elde).(co).(am).(itf).eta;
        P(irev, i) = sum(n * F .* vols .* Rvol .* eta);

        if model.(elde).(co).(am).(itf).includeEntropyChange
            dUdT = st.(elde).(co).(am).(itf).dUdT;
            P(irevb, i) = sum(n * F .* vols .* Rvol .* T(cmap) .* dUdT);
        end
    end

    % Diffusion (electrolyte)
    P(10, i) = sum(st.ThermalModel.jHeatChemicalSource);
end

Pres = Ptot - sum(P, 1);
time_s = reshape(time, 1, []);
t_h = time_s / hour;

% Accumulated energy (J)
Eparts = zeros(size(P));
for k = 1 : nlabels
    Eparts(k, :) = cumtrapz(time_s, P(k, :));
end
Eres = cumtrapz(time_s, Pres);
Etot = cumtrapz(time_s, Ptot);

figure(4); clf;
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% Match Julia plotting structure: stacked signed areas with residual included.
labels_plot = [labels, {'Residual'}];
colors_plot = [colors; 0.45 0.45 0.45];
Pplot = [P; Pres];
Eplot = [Eparts; Eres];
nplot = numel(labels_plot);

% Top: total heat production + stacked signed areas
nexttile; hold on; box on;
hbands = gobjects(nplot, 1);
cum_pos = zeros(1, nt);
cum_neg = zeros(1, nt);
for k = 1 : nplot
    v = Pplot(k, :);
    y1 = zeros(1, nt);
    y2 = zeros(1, nt);
    ip = (v >= 0);
    in = ~ip;
    y1(ip) = cum_pos(ip);
    y2(ip) = cum_pos(ip) + v(ip);
    cum_pos(ip) = y2(ip);
    y1(in) = cum_neg(in);
    y2(in) = cum_neg(in) + v(in);
    cum_neg(in) = y2(in);
    xv = [t_h, fliplr(t_h)];
    yv = [y1, fliplr(y2)];
    hbands(k) = patch(xv, yv, colors_plot(k, :), ...
                      'EdgeColor', 'none', ...
                      'FaceAlpha', 0.65, ...
                      'DisplayName', labels_plot{k});
end
h_total = plot(t_h, Ptot, 'k--', 'LineWidth', 3, 'DisplayName', 'Total');
title('Total heat production', 'FontSize', 30, 'FontWeight', 'bold');
xlabel('Time  /  h', 'FontSize', 24);
ylabel('Total heat production  /  W', 'FontSize', 24);
set(gca, 'FontSize', 18);
legend([hbands; h_total], ...
       [labels_plot, {'Total'}], ...
       'Location', 'eastoutside', 'FontSize', 16);

% Bottom: accumulated source contributions + stacked signed areas
nexttile; hold on; box on;
hbands2 = gobjects(nplot, 1);
cum_pos = zeros(1, nt);
cum_neg = zeros(1, nt);
for k = 1 : nplot
    v = Eplot(k, :);
    y1 = zeros(1, nt);
    y2 = zeros(1, nt);
    ip = (v >= 0);
    in = ~ip;
    y1(ip) = cum_pos(ip);
    y2(ip) = cum_pos(ip) + v(ip);
    cum_pos(ip) = y2(ip);
    y1(in) = cum_neg(in);
    y2(in) = cum_neg(in) + v(in);
    cum_neg(in) = y2(in);
    xv = [t_h, fliplr(t_h)];
    yv = [y1, fliplr(y2)];
    hbands2(k) = patch(xv, yv, colors_plot(k, :), ...
                       'EdgeColor', 'none', ...
                       'FaceAlpha', 0.65, ...
                       'DisplayName', labels_plot{k});
end
h_total_acc = plot(t_h, Etot, 'k--', 'LineWidth', 3, 'DisplayName', 'Accumulated/total');
title('Accumulation of individual heat sources', 'FontSize', 30, 'FontWeight', 'bold');
xlabel('Time  /  h', 'FontSize', 24);
ylabel('Accumulated heat  /  J', 'FontSize', 24);
set(gca, 'FontSize', 18);
legend([hbands2; h_total_acc], ...
       [labels_plot, {'Accumulated/total'}], ...
       'Location', 'eastoutside', 'FontSize', 16);


% ============================================================
% Store heat source labels with their (time, position) matrices
% ============================================================

customKeys = { ...
    "OhmicCurrentCollector (pos)", ...
    "OhmicCurrentCollector (neg)", ...
    "OhmicActiveMaterial (pos)", ...
    "OhmicActiveMaterial (neg)", ...
    "OhmicElectrolyte (elyte)", ...
    "ReactionReversible (pos)", ...
    "ReactionReversible (neg)", ...
    "ReactionIrreversible (pos)", ...
    "ReactionIrreversible (neg)", ...
    "DiffusionElectrolyte (elyte)" ...
};

% Create a struct with valid MATLAB fieldnames
heatSourceDataStruct = struct();

for k = 1:nlabels
    key = customKeys{k};
    field = matlab.lang.makeValidName(key);  % convert to valid struct name
    heatSourceDataStruct.(field) = [time_s(:), P(k, :).'];
end

% Add Residual and Total for completeness
heatSourceDataStruct.Residual = [time_s(:), Pres(:)];
heatSourceDataStruct.Total    = [time_s(:), Ptot(:)];

assignin('base', 'heatSourceData', heatSourceDataStruct);


%% Export comparison bundle for BattMo.jl
% Save variables needed by Julia-side thermal parity checks, including
% effective thermal vectors from MATLAB thermal model setup.
mat_out_candidates = { ...
    fullfile('test', 'data', 'matlab_files', 'run_only_thermal.mat'), ...
    fullfile('..', 'BattMo.jl', 'test', 'data', 'matlab_files', 'run_only_thermal.mat'), ...
    fullfile('..', '..', 'BattMo.jl', 'test', 'data', 'matlab_files', 'run_only_thermal.mat'), ...
    fullfile('..', '..', '..', 'BattMo.jl', 'test', 'data', 'matlab_files', 'run_only_thermal.mat') ...
};

output_isothermal_states = output_isothermal.states;
output_isothermal_time = output_isothermal.time;

vars_to_save = {'time', 'E', 'sourceTerms', 'states_thermal', ...
                'output_isothermal_states', 'output_isothermal_time', ...
                'effectiveThermalConductivity', 'effectiveVolumetricHeatCapacity', ...
                'thermalBoundaryNeighbors', 'thermalBoundaryAreas', 'jsonstruct', 'inputparams_thermal', 'states_heat','heatSourceData'};

saved_any = false;
for iout = 1:numel(mat_out_candidates)
    mat_out = mat_out_candidates{iout};
    out_dir = fileparts(mat_out);
    if exist(out_dir, 'dir')
        save(mat_out, vars_to_save{:}, '-v7.3');
        fprintf('Saved MATLAB thermal comparison data to: %s\n', mat_out);
        saved_any = true;
    end
end

if ~saved_any
    warning('Could not find a target directory to save run_only_thermal.mat');
end



