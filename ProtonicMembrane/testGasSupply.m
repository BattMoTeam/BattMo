clear all

mrstModule add ad-core

filename = '/home/xavier/Matlab/Projects/battmo/ProtonicMembrane/gaslayersupply.json';
jsonstruct = fileread(filename);
jsonstruct = jsondecode(jsonstruct);

paramobj = ProtonicMembraneGasSupplyInputParams(jsonstruct);

nx = 10;
ny = 7;

G = cartGrid([nx, ny]);
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

T = 10;
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

