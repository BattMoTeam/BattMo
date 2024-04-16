clear all

mrstModule add ad-core mrst-gui

filename = '/home/xavier/Matlab/Projects/battmo/ProtonicMembrane/gas_supply.json';
jsonstruct = fileread(filename);
jsonstruct = jsondecode(jsonstruct);

%%
Dmult = 1e-7;
fprintf('Diffusion coefficient multiplier = %g\n', Dmult);
jsonstruct.diffusionCoefficients = Dmult*jsonstruct.diffusionCoefficients;

%%
rate = 1e-3;
fprintf('Use rate = %g\n', rate);
jsonstruct.control(1).values(1) = rate;

%%
permMult = 1e-3;
fprintf('Use permMult = %g\n', permMult);
jsonstruct.permeability = permMult*jsonstruct.permeability;

inputparams = ProtonicMembraneGasSupplyInputParams(jsonstruct);

dimcase = '1D';

switch dimcase

  case '2D'

    gen = GasSupplyGridGenerator2D();

    gen.nx = 10;
    gen.ny = 1000;
    gen.lx = 0.5*milli*meter;
    gen.ly = 1.5*milli*meter;
    
  case '1D'
    
    gen = GasSupplyGridGenerator1D();

    gen.nx = 1000;
    gen.lx = 1.5*milli*meter;
    gen.faceArea = 0.5*milli*meter;
    
  otherwise

    error('dimcase not recognized')
end

inputparams = gen.updateInputParams(inputparams);

% Setup model

model = ProtonicMembraneGasSupply(inputparams);
model = model.setupForSimulation();

cgt = model.cgt;
cgp = model.cgp;


%% Setup initial state

initstate = model.setupInitialState();

%%
mfInit = 0.2;
fprintf('Use mfInit = %g\n', mfInit);
initstate.massfractions{1}             = mfInit + 0*initstate.massfractions{1};
initstate.GasSupplyBc.massfractions{1} = mfInit + 0*initstate.GasSupplyBc.massfractions{1};

%% Setup scalings

gasInd = model.gasInd;

pH2O = initstate.pressure(1);
rho  = initstate.density(1);

%% pick-up rate value in control
%
comptypecouptbl = model.helpers.comptypecouptbl;
ctrlvals        = model.helpers.ctrlvals;

clear comptypetbl2
comptypetbl2.type = [2]; % type=2 for rate control
comptypetbl2.comp = [1]; % component index = 1 for rate value
comptypetbl2 = IndexArray(comptypetbl2);

map = TensorMap();
map.fromTbl = comptypecouptbl;
map.toTbl = comptypetbl2;
map.mergefds = {'comp', 'type'};
map = map.setup();

rate = max(map.eval(ctrlvals));

scalFlux     = rate/gen.nx;
scalPressure = pH2O;

model.scalings = {{{'massConses', 1}, scalFlux}, ...
                  {{'massConses', 2}, scalFlux}, ...
                  {{'Control', 'pressureEq'}, scalPressure}, ...
                  {{'Control', 'rateEq'}, gen.nx*scalFlux}, ...
                  {{'GasSupplyBc', 'bcFluxEquations', 1}, scalFlux}, ...
                  {{'GasSupplyBc', 'bcFluxEquations', 2}, scalFlux}};

T = 1e-3*second;
N = 1;
dt = rampupTimesteps(T, T/N, 1, 'threshold_error', 1e-15);

step.val = dt;
step.control = ones(numel(step.val), 1);

control.src = [];

schedule = struct('control', control, 'step', step);

nls = NonLinearSolver();
nls.maxIterations = 20;
nls.errorOnFailure = false;
nls.verbose = true;

model.verbose = true;

model.nonlinearTolerance = 1e-5;

[~, states, report] = simulateScheduleAD(initstate, model, schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls);

%%

figure
plotToolbar(model.grid, states);
% caxis([0.2, 0.4])
% uit = findobj(gcf, 'Tooltip', 'Freeze caxis');
% uit.State = 'on';



