filename = fullfile(battmoDir(), 'Diffusion', 'Parameters', 'gasdiffusioncell.json');
jsonstruct = parseBattmoJson(filename);
inputparams = GasDiffusionCellInputParams(jsonstruct);

gen = GasDiffusionCellGridGenerator1D(1, 10);
inputparams = gen.updateGasDiffusionCellInputParams(inputparams);

model = GasDiffusionCell(inputparams);
model = model.equipModelForComputation();

initstate = model.setupInitialState(jsonstruct);

cgt = model.cgt;
cgt.printRootVariables;
cgt.printTailVariables;

% cgp = model.cgp;

%% Setup initial state

initstate = model.setupInitialState(jsonstruct);

%% Setup control
nsteps   = 10;
stepsize = 1*second;

clear step
step.val     = stepsize*ones(nsteps, 1);
step.control = ones(nsteps, 1);
control.src = [];

schedule = struct('control', control, 'step', step);

nls = NonLinearSolver();
nls.errorOnFailure = false;

model.nonlinearTolerance = 1e-2;

%% Run simulation

model.verbose = true;
[~, states, report] = simulateScheduleAD(initstate, model, schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls);


