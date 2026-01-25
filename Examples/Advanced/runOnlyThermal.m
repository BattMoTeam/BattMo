%% Decoupled electro-chemical and thermal simulation

%% setup material property input
jsonfilename = fullfile('ParameterData'        , ...
                        'BatteryCellParameters', ...
                        'LithiumIonBatteryCell', ...
                        'lithium_ion_battery_nmc_graphite.json');
jsonstruct_material = parseBattmoJson(jsonfilename);

jsonstruct_material.include_current_collectors = true;

%% Setup geometry input
%
% We use a simple 3d-geometry
jsonfilename = fullfile('Examples'     , ...
                        'JsonDataFiles', ...
                        'geometry3d.json');
jsonstruct_geometry = parseBattmoJson(jsonfilename);

%% Setup  Control input
%
jsonfilename = fullfile('Examples', 'JsonDataFiles', 'cc_discharge_control.json');
jsonstruct_control = parseBattmoJson(jsonfilename);

%% Setup full input
%
jsonstruct = mergeJsonStructs({jsonstruct_geometry , ...
                               jsonstruct_material , ...
                               jsonstruct_control}, 'warn', false);


jsonstruct.Electrolyte.ionicConductivity.functionName = 'modified_electrolyte_conductivity';

%% Run non thermal simulation

jsonstruct_isothermal = setJsonStructField(jsonstruct, 'use_thermal', false, 'handleMisMatch', 'quiet');
output_isothermal = runBatteryJson(jsonstruct_isothermal);

states = output_isothermal.states;
times  = output_isothermal.time;

%% Setup model with thermal support

output_fullycoupled = runBatteryJson(jsonstruct);
model = output_fullycoupled.model;

% We recover the initial state from which we will extract the temperature (it is actually a constant value...)
initstate = model.setupInitialState(jsonstruct);

%% Add source terms to the state output using the thermal model

for istate = 1 : numel(states)
    states{istate} = model.evalVarName(states{istate}, {'ThermalModel', 'jHeatSource'});
end

%% Setup source term helper object
% We use the thermal source values in states that we just computted

sourceTerms = cellfun(@(state) state.ThermalModel.jHeatSource, states, 'uniformoutput', false);

hss = HeatSourceSetup(sourceTerms, times);

%% Setup thermal model only

inputparams_thermal = output_fullycoupled.inputparams.ThermalModel;

inputparams_thermal.effectiveThermalConductivity    = output_fullycoupled.model.ThermalModel.effectiveThermalConductivity;
inputparams_thermal.effectiveVolumetricHeatCapacity = output_fullycoupled.model.ThermalModel.effectiveVolumetricHeatCapacity;

model_thermal = ThermalComponent(inputparams_thermal);
model_thermal.isRootSimulationModel = true;
model_thermal = model_thermal.equipModelForComputation();

%% Setup the schedule

times = hss.times;

clear step
step.val     = [times(1); diff(times)];
step.control = ones(numel(times), 1);

clear control
control.src = @(time) hss.eval(time);

schedule = struct('control', control, ...
                  'step'   , step);

clear state0;
state0.T = initstate.ThermalModel.T;

simInput = struct('model'    , model_thermal, ...
                  'initstate', state0, ...
                  'schedule' , schedule);

simsetup = SimulationSetup(simInput);

states_thermal = simsetup.run();

temp_plot_thermal

