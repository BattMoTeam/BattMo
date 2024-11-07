
dosetup = false;

if dosetup
    states      = output.states;
    model       = output.model;
    inputparams = output.inputparams;
    times       = output.time;
    for istate = 1 : numel(states)
        states{istate} = model.addVariables(states{istate});
    end
    initstate = model.setupInitialState();

    inputparams = inputparams.ThermalModel;

end

sourceTerms = cellfun(@(state) state.ThermalModel.jHeatSource, states, 'uniformoutput', false);

hss = HeatSourceSetup(sourceTerms, times);

inputparams.effectiveThermalConductivity    = output.model.ThermalModel.effectiveThermalConductivity;
inputparams.effectiveVolumetricHeatCapacity = output.model.ThermalModel.effectiveVolumetricHeatCapacity;

model = ThermalComponent(inputparams);
model.isRootSimulationModel = true;
model = model.equipModelForComputation();

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

model.verbose = true;
[~, states, report] = simulateScheduleAD(state0, model, schedule);


