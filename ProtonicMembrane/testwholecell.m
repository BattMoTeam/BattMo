clear all

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

jsonstruct.(gs).control(1).values = [pH2O; pO2];
jsonstruct.(gs).control(2).values = [pH2O; pO2];

paramobj = ProtonicMembraneCellWithGasSupplyInputParams(jsonstruct);

gen = GasSupplyPEMgridGenerator2D();

gen.nxCell      = 10;
gen.nxGasSupply = 10;
gen.lxCell      = 22*micro*meter;
gen.lxGasSupply = 1*milli*meter;

gen.ny = 10;
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

%% Setup schedule

tswitch = 1;
T       = 2; % This is not a real time scale, as all the model deals with equilibrium

N1  = 20;
dt1 = tswitch/N1;
N2  = 20;
dt2 = (T - tswitch)/N2;

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

model.nonlinearTolerance = 1e-8;
model.verbose = true;

%% Start simulation

[~, states, report] = simulateScheduleAD(initstate, model, schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls); 




