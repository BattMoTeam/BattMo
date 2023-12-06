clear all

mrstDebug(20);

mrstModule add ad-core mrst-gui

gs    = 'GasSupply';
ce    = 'Cell';
an    = 'Anode';
ct    = 'Cathode';
elyte = 'Electrolyte';
ctrl  = 'Control';
            
clear jsonstruct

filename = 'ProtonicMembrane/gas_supply.json';
jsonstruct.GasSupply = parseBattmoJson(filename);
filename = 'ProtonicMembrane/protonicMembrane.json';
jsonstruct.Cell = parseBattmoJson(filename);

pH2O = 475000;
pO2  = 112500;

p = pH2O + pO2;
s = pH2O/p;

jsonstruct.(gs).control(1).values = [1.1*p; s];
jsonstruct.(gs).control(2).values = [p; s];

paramobj = ProtonicMembraneCellWithGasSupplyInputParams(jsonstruct);

gen = GasSupplyPEMgridGenerator2D();

gen.nxCell      = 1000;
gen.nxGasSupply = 50;
gen.lxCell      = 22*micro*meter;
gen.lxGasSupply = 1*milli*meter;

gen.ny = 50;
gen.ly = 1;

paramobj = gen.updateInputParams(paramobj);

doplot = false;

if doplot
    
    close all

    figure('position', [337, 757, 3068, 557])
    plotGrid(paramobj.G)
    plotGrid(paramobj.Cell.G, 'facecolor', 'red')
    plotGrid(paramobj.GasSupply.G, 'facecolor', 'blue')

    plotGrid(paramobj.GasSupply.G, paramobj.GasSupply.couplingTerms{1}.couplingcells);
    plotGrid(paramobj.GasSupply.G, paramobj.GasSupply.couplingTerms{2}.couplingcells);
    plotGrid(paramobj.Cell.G, paramobj.Cell.couplingTerms{1}.couplingcells(:, 2));
    plotGrid(paramobj.Cell.G, paramobj.Cell.couplingTerms{2}.couplingcells(:, 2));
    plotGrid(paramobj.GasSupply.G, paramobj.GasSupply.couplingTerms{3}.couplingcells );

    return
    
end

model = ProtonicMembraneCellWithGasSupply(paramobj);

model = model.setupForSimulation();

%% Setup initial state

initstate = model.setupInitialState();

%% setup scalings

gasInd = model.(gs).gasInd;

pH2O         = initstate.(gs).pressure(1);
rho          = initstate.(gs).density(1);
scalFlux     = 1/gen.nxGasSupply*rho*model.(gs).permeability/model.(gs).viscosity*pH2O/gen.ly;
scalPressure = pH2O;

model.scalings = {{{gs, 'massConses', 1}, scalFlux}, ...
                  {{gs, 'massConses', 2}, scalFlux}, ...
                  {{gs, 'GasSupplyBc', 'controlEquations'}, scalPressure}, ...
                  {{gs, 'GasSupplyBc', 'boundaryEquations', 1}, scalFlux}, ...
                  {{gs, 'GasSupplyBc', 'boundaryEquations', 2}, scalFlux}};

drivingforces.src = @(time) controlfunc(time, 0, 1, 2, 'order', 'I-first');
initstate = model.evalVarName(initstate, {ce, elyte, 'sigmaEl'}, {{'drivingForces', drivingforces}});
initstate = model.evalVarName(initstate, {ce, elyte, 'sigmaHp'}, {{'drivingForces', drivingforces}});

sigmaHp = initstate.(ce).(elyte).sigmaHp(1);
sigmaEl = initstate.(ce).(elyte).sigmaEl(1);
phi0    = abs(model.(ce).(an).E_0 - model.(ce).(ct).E_0); % characteristic voltage
T       = model.(ce).(elyte).operators.T_all(1);

sHp = T*sigmaHp*phi0;
sEl = T*sigmaEl*phi0;

model.scalings =  horzcat(model.scalings, ...
                          {{{ce, elyte, 'massConsHp'}     , sHp}, ...
                           {{ce, elyte, 'chargeConsEl'}   , sEl}, ...
                           {{ce, an   , 'chargeCons'}     , sEl}, ...
                           {{ce, an   , 'iElEquation'}    , sEl}, ...
                           {{ce, an   , 'iHpEquation'}    , sHp}, ...
                           {{ce, ct   , 'chargeCons'}     , sEl}, ...
                           {{ce, ct   , 'iElEquation'}    , sEl}, ...
                           {{ce, ct   , 'iHpEquation'}    , sHp}, ...
                           {{ce, ctrl , 'controlEquation'}, sEl}, ...
                           {{ce, 'anodeChargeCons'}       , sEl}});

%% Setup schedule

tswitch = 0.5;
T       = 2*hour; 

N1  = 20;
dt1 = tswitch*T/N1;
N2  = 20;
dt2 = T*(1 - tswitch)/N2;

step.val = [dt1*ones(N1, 1); dt2*ones(N2, 1)];
step.control = ones(numel(step.val), 1);

Imax = 0;

control.src = @(time) controlfunc(time, Imax, tswitch, T, 'order', 'I-first');

schedule = struct('control', control, 'step', step); 

%% Setup nonlinear solver

nls                = NonLinearSolver();
nls.maxIterations  = 20;
nls.errorOnFailure = false;
nls.verbose        = true;

model.nonlinearTolerance = 1e-5;
model.verbose = true;

%% Start simulation

dopack = true;
clearSimulation = false;

if dopack

    name = 'testwholecell';
    problem = packSimulationProblem(initstate, model, schedule, [], ...
                                    'Name'           , name      , ...
                                    'NonLinearSolver', nls);
    if clearSimulation
        %% clear previously computed simulation
        clearPackedSimulatorOutput(problem, 'prompt', false);
    end
    simulatePackedProblem(problem);
    [globvars, states, reports] = getPackedSimulatorOutput(problem);
    
else
    
    [~, states, report] = simulateScheduleAD(initstate, model, schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls); 

end

%%

close all

set(0, 'defaultlinelinewidth', 3);
set(0, 'defaultaxesfontsize', 15);

N = gen.nxCell;
xc = model.(ce).(elyte).G.cells.centroids(1 : N, 1);

state = states{end};

X = reshape(model.(ce).(elyte).G.cells.centroids(:, 1), N, []);
Y = reshape(model.(ce).(elyte).G.cells.centroids(:, 2), N, []);

figure
val = state.(ce).(elyte).pi;
Z = reshape(val, N, []);
surf(X, Y, Z, 'edgecolor', 'none');
title('pi')
xlabel('x [m]')
view(20, 15)
colorbar

figure
val = state.(ce).(elyte).pi - state.(ce).(elyte).phi;
Z = reshape(val, N, []);
surf(X, Y, Z, 'edgecolor', 'none');
title('E')
xlabel('x [m]')
view(20, 15)
colorbar

figure
val = state.(ce).(elyte).phi;
Z = reshape(val, N, []);
surf(X, Y, Z, 'edgecolor', 'none');
title('phi')
xlabel('x [m]')
view(20, 15)
colorbar




