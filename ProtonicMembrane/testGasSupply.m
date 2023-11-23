clear all

mrstModule add ad-core mrst-gui

filename = '/home/xavier/Matlab/Projects/battmo/ProtonicMembrane/gas_supply.json';
jsonstruct = fileread(filename);
jsonstruct = jsondecode(jsonstruct);

paramobj = ProtonicMembraneGasSupplyInputParams(jsonstruct);
gen = GasSupplyGridGenerator2D();

gen.nx = 100;
gen.ny = 70;
gen.lx = 10;
gen.ly = 10;

paramobj = gen.updateInputParams(paramobj);

% Setup model

model = ProtonicMembraneGasSupply(paramobj);
model = model.setupStandAlone();

cgt = model.computationalGraph;

% Setup initial state

nc  = model.G.cells.num;
nbc = model.GasSupplyBc.getNumberBcFaces();

pH2O = jsonstruct.control(2).values(1);
pO2  = jsonstruct.control(2).values(2);

gasInd = model.gasInd;

initstate.pressures{gasInd.H2O}              = pH2O*ones(nc, 1);
initstate.pressures{gasInd.O2}               = pO2 *ones(nc, 1);
initstate.GasSupplyBc.pressures{gasInd.H2O}  = pH2O*ones(nbc, 1);
initstate.GasSupplyBc.pressures{gasInd.O2}   = pO2 *ones(nbc, 1);
initstate.GasSupplyBc.massFluxes{gasInd.H2O} = zeros(nbc, 1);
initstate.GasSupplyBc.massFluxes{gasInd.O2}  = zeros(nbc, 1);

initstate = model.evalVarName(initstate, VarName({}, 'densities', model.nGas));

%% setup scalings

rho = initstate.densities{1}(1);
scalFlux     = 1/gen.nx*rho*model.permeability/model.viscosity*pH2O/gen.ly;
scalPressure = pH2O;

model.scalings = {{{'massConses', 1}, scalFlux}, ...
                  {{'massConses', 2}, scalFlux}, ...
                  {{'GasSupplyBc', 'controlEquation'}, scalPressure}, ...
                  {{'GasSupplyBc', 'boundaryEquations', 1}, scalFlux}, ...
                  {{'GasSupplyBc', 'boundaryEquations', 2}, scalFlux}};

T = 1*minute;
N = 10;

dt = T/N;

step.val = dt*ones(N, 1);
step.control = ones(numel(step.val), 1);

control.src = [];

schedule = struct('control', control, 'step', step);

nls = NonLinearSolver();
nls.maxIterations = 20;
nls.errorOnFailure = false;
nls.verbose = true;

model.nonlinearTolerance = 1e-8;

[~, states, report] = simulateScheduleAD(initstate, model, schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls);

%%

close all

ind = cellfun(@(state) ~isempty(state), states);
states = states(ind);


states = cellfun(@(state) model.evalVarName(state, 'pressure'), states, 'un', false);
plotToolbar(model.G, states);
