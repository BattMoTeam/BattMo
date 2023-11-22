clear all

mrstModule add ad-core mrst-gui

filename = '/home/xavier/Matlab/Projects/battmo/ProtonicMembrane/gaslayersupply.json';
jsonstruct = fileread(filename);
jsonstruct = jsondecode(jsonstruct);

paramobj = ProtonicMembraneGasSupplyInputParams(jsonstruct);

nx = 100;
ny = 70;
lx = 10;
ly = 10;

G = cartGrid([nx, ny], [lx, ly]);
G = computeGeometry(G);

tbls = setupSimpleTables(G);
cellfacetbl = tbls.cellfacetbl;

clear bcfacecouptbl1;
bcfacecouptbl1.faces = (nx + 1)*ny + (1 : floor(nx/2))';
bcfacecouptbl1 = IndexArray(bcfacecouptbl1);
nc = bcfacecouptbl1.num;
bcfacecouptbl1 =  bcfacecouptbl1.addInd('coup', ones(nc, 1));

clear bcfacecouptbl2;
bcfacecouptbl2.faces = (nx + 1)*ny + ny*nx + (floor(nx/2) + 1 : nx)';
bcfacecouptbl2 = IndexArray(bcfacecouptbl2);
nc = bcfacecouptbl2.num;
bcfacecouptbl2 =  bcfacecouptbl2.addInd('coup', 2*ones(nc, 1));

bcfacecouptbl = concatIndexArray(bcfacecouptbl1, bcfacecouptbl2, {});

bccellfacecouptbl = crossIndexArray(bcfacecouptbl, cellfacetbl, {'faces'});

coupnames = {'input', 'output'};

couplingTerms = {};

for  icoup = 1 : numel(coupnames)

    coupname = coupnames{icoup};
    name = sprintf('External %s', coupname);
    compnames = sprintf('External-%s', coupname);
    coupTerm = couplingTerm(name, compnames);

    clear couptbl;
    couptbl.coup =  icoup;
    couptbl = IndexArray(couptbl);

    tbl =  crossIndexArray(couptbl, bccellfacecouptbl, {'coup'});

    coupTerm.couplingfaces = tbl.get('faces');
    coupTerm.couplingcells = tbl.get('cells');

    couplingTerms{end + 1} = coupTerm;
    
end

paramobj.couplingTerms = couplingTerms;
paramobj.G = G;

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
scalFlux     = 1/nx*rho*model.permeability/model.viscosity*pH2O/ly;
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
