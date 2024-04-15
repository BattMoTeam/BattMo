clear all

mrstModule add ad-core mrst-gui

filename = '/home/xavier/Matlab/Projects/battmo/ProtonicMembrane/gas_supply.json';
jsonstruct = fileread(filename);
jsonstruct = jsondecode(jsonstruct);

%%
Dmult = 1;
fprintf('Diffusion coefficient multiplier = %g\n', Dmult);
jsonstruct.diffusionCoefficients = Dmult*jsonstruct.diffusionCoefficients;

%%
rate = 1e-4;
fprintf('Use rate = %g\n', rate);
jsonstruct.control(1).values(1) = rate;

%%
perm = 0;
fprintf('Use perm = %g\n', perm);
jsonstruct.permeability = perm;

inputparams = ProtonicMembraneGasSupplyInputParams(jsonstruct);

gen = GasSupplyGridGenerator2D();

gen.nx = 100;
gen.ny = 70;
gen.lx = 0.5*milli*meter;
gen.ly = 1.5*milli*meter;

inputparams = gen.updateInputParams(inputparams);

% Setup model

model = ProtonicMembraneGasSupply(inputparams);
model = model.setupForSimulation();

cgt = model.cgt;
cgp = model.cgp;


%% Setup initial state

initstate = model.setupInitialState();

%% Setup scalings

gasInd = model.gasInd;

pH2O         = initstate.pressure(1);
rho          = initstate.density(1);
scalFlux     = 1/gen.nx*(rho*model.permeability/model.viscosity*pH2O + rho*model.diffusionCoefficients(1))/gen.ly;
scalPressure = pH2O;

model.scalings = {{{'massConses', 1}, scalFlux}, ...
                  {{'massConses', 2}, scalFlux}, ...
                  {{'Control', 'pressureEq'}, scalPressure}, ...
                  {{'Control', 'rateEq'}, gen.nx*scalFlux}, ...
                  {{'GasSupplyBc', 'bcFluxEquations', 1}, scalFlux}, ...
                  {{'GasSupplyBc', 'bcFluxEquations', 2}, scalFlux}};


T = 1*second;
N = 1;
dt = rampupTimesteps(T, T/N, 10, 'threshold_error', 1e-15);

step.val = dt;
step.control = ones(numel(step.val), 1);

control.src = [];

schedule = struct('control', control, 'step', step);

nls = NonLinearSolver();
nls.maxIterations = 20;
nls.errorOnFailure = false;
nls.verbose = true;

model.verbose = true;

model.nonlinearTolerance = 1e-8;

[~, states, report] = simulateScheduleAD(initstate, model, schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls);

%%

figure
plotToolbar(model.grid, states);
caxis([0.2, 0.4])
% uit = findobj(gcf, 'Tooltip', 'Freeze caxis');
% uit.State = 'on';



