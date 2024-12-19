clear all

% mrstDebug(20);

mrstModule add ad-core mrst-gui

gs    = 'GasSupply';
ce    = 'Electrolyser';
an    = 'Anode';
ct    = 'Cathode';
elyte = 'Electrolyte';
ctrl  = 'Control';
            
clear jsonstruct_material

filename = 'ProtonicMembrane/jsonfiles/gas-supply-whole-cell.json';
jsonstruct_material.GasSupply = parseBattmoJson(filename);
filename = 'ProtonicMembrane/jsonfiles/protonicMembrane.json';
jsonstruct_material.Electrolyser = parseBattmoJson(filename);

filename = 'ProtonicMembrane/jsonfiles/2d-cell-geometry.json';
jsonstruct_geometry = parseBattmoJson(filename);

jsonstruct_material = removeJsonStructFields(jsonstruct_material, ...
                                             {'Electrolyser', 'Electrolyte', 'Nx'}, ...
                                             {'Electrolyser', 'Electrolyte', 'xlength'});

jsonstruct = mergeJsonStructs({jsonstruct_material, ...
                               jsonstruct_geometry});

%% Adjust diffusion
Dmult = 1e-3;
fprintf('Diffusion coefficient multiplier : %g\n', Dmult);
jsonstruct.(gs).diffusionCoefficients = Dmult*jsonstruct.(gs).diffusionCoefficients;
%%

%% Adjust rate
rate = 1e-6;
fprintf('Rate value : %g\n', rate);
jsonstruct.(gs).control(1).values(1) = rate;

%%

inputparams        = ProtonicMembraneCellInputParams(jsonstruct);
[inputparams, gen] = setupProtonicMembraneCellGrid(inputparams, jsonstruct);

%% Adjust I

I = 0.5/((centi*meter)^2) * gen.ly;
fprintf('use I = %g\n', I);

%%

rungaslayer = false;

if rungaslayer

    rungaslayer_only;
    
end

doinitialisation = false;

if doinitialisation

    run_initialization;
    
else

    switch I
      case 0
        molfluxref = 1.7236e-06; % value for I = 0
      case 7.5
        molfluxref = 7.35713e-05; % mol/second, for I = 7.5
      otherwise
        error('value of molflux ref did not appear to have been computed precedently');
    end
    
end

%%

model = ProtonicMembraneCell(inputparams);

figure('position', [337, 757, 3068, 557])
plotGrid(model.grid)
plotGrid(model.Electrolyser.grid, 'facecolor', 'red')
plotGrid(model.GasSupply.grid, 'facecolor', 'blue')

plotGrid(model.GasSupply.grid, model.GasSupply.couplingTerms{1}.couplingcells);
plotGrid(model.GasSupply.grid, model.GasSupply.couplingTerms{2}.couplingcells);
plotGrid(model.Electrolyser.grid, model.Electrolyser.couplingTerms{1}.couplingcells(:, 2));
plotGrid(model.Electrolyser.grid, model.Electrolyser.couplingTerms{2}.couplingcells(:, 2));
plotGrid(model.GasSupply.grid, model.couplingTerm.couplingcells(:, 1) );

model = model.setupForSimulation();


%% Print cutoff values
%

fprintf('Cutoff min pressure : %g\n'     , model.pmin);
fprintf('Cutoff min mass fraction : %g\n', model.mfmin);
fprintf('Cutoff max mass fraction : %g\n', model.mfmax);
fprintf('Cutoff abs potentials : %g\n'   , model.phimax);

%% Setup initial state

initstate = model.setupInitialState();

%% Setup scalings

gasInd = model.(gs).gasInd;

pH2O         = initstate.(gs).pressure(1);
rho          = initstate.(gs).density(1);
scalFlux     = molfluxref/model.(gs).molecularWeights(1)/gen.ny; % in kg/s
scalPressure = pH2O;

model.scalings = {{{gs, 'massConses', 1}, scalFlux}                      , ...
                  {{gs, 'massConses', 2}, scalFlux}                      , ...
                  {{gs, 'Control', 'pressureEq'}, scalPressure}          , ...
                  {{gs, 'Control', 'rateEq'}, gen.nxGasSupply*scalFlux}           , ...
                  {{gs, 'GasSupplyBc', 'bcFluxEquations', 1}, scalFlux}, ...
                  {{gs, 'GasSupplyBc', 'bcFluxEquations', 2}, scalFlux}};

drivingforces.src = @(time) controlfunc(time, 0, 1, 2, 'order', 'I-first');
initstate = model.evalVarName(initstate, {ce, elyte, 'sigmaEl'}, {{'drivingForces', drivingforces}});
initstate = model.evalVarName(initstate, {ce, elyte, 'sigmaHp'}, {{'drivingForces', drivingforces}});

sigmaHp = initstate.(ce).(elyte).sigmaHp(1);
sigmaEl = initstate.(ce).(elyte).sigmaEl(1);
phi0    = abs(model.(ce).(an).E_0 - model.(ce).(ct).E_0); % characteristic voltage
T       = model.(ce).(elyte).G.getTrans();
T       = T(1);

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

tswitch   = 0.1;
totaltime = 10*minute;

timeswitch = tswitch*totaltime;

%% adjust N2
N1  = 10;
fprintf('use N1 = %g\n', N1);
%%

dt1 = timeswitch/N1;
t1  = linspace(0, dt1, N1 + 1)';
alpha = 40;
t1 = log(alpha*t1 + 1)./log(alpha*dt1 + 1)*dt1;

steps1 = diff(t1);

%% adjust N2
N2  = 20;
fprintf('use N2 = %g\n', N2);
%%

dt2 = (totaltime - timeswitch)/N2;
steps2 = rampupTimesteps(totaltime -timeswitch, dt2, 5);

%%

step.val = [steps1; steps2];
step.control = ones(numel(step.val), 1);

control.src = @(time) controlfunc(time, I, timeswitch, totaltime, 'order', 'alpha-equal-beta');

schedule = struct('control', control, 'step', step); 

%% Setup nonlinear solver

ts = IterationCountTimeStepSelector('targetIterationCount', 5, ...
                                    'verbose', true);

nls = NonLinearSolver();
nls.timeStepSelector = ts;
nls.maxIterations    = 15;
nls.errorOnFailure   = false;
nls.verbose          = true;

model.nonlinearTolerance = 1e-5;
model.verbose = true;

%% Start simulation

dopack = true;
clearSimulation = true;

if dopack

    name = 'testwholecell';
    problem = packSimulationProblem(initstate, model, schedule, []             , ...
                                    'ExtraArguments', {'OutputMinisteps', true}, ...
                                    'Name'           , name                    , ...
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

%% plotting

close all

set(0, 'defaultlinelinewidth', 3);
set(0, 'defaultaxesfontsize', 15);

N = gen.nxElectrolyser;
xc = model.(ce).(elyte).grid.cells.centroids(1 : N, 1);

state = states{end};

state = model.addVariables(state);

X = reshape(model.(ce).(elyte).grid.cells.centroids(:, 1), N, [])/(milli*meter);
Y = reshape(model.(ce).(elyte).grid.cells.centroids(:, 2), N, [])/(milli*meter);

figure
val = state.(ce).(elyte).pi;
Z = reshape(val, N, []);
surf(X, Y, Z, 'edgecolor', 'none');
title('pi / V')
xlabel('x / mm')
ylabel('y / mm')
view(45, 31)
colorbar

figure
val = state.(ce).(elyte).pi - state.(ce).(elyte).phi;
Z = reshape(val, N, []);
surf(X, Y, Z, 'edgecolor', 'none');
title('E / V')
xlabel('x / mm')
ylabel('y / mm')
view(45, 31)
colorbar

figure
val = state.(ce).(elyte).phi;
Z = reshape(val, N, []);
surf(X, Y, Z, 'edgecolor', 'none');
title('phi / V')
xlabel('x / mm')
ylabel('y / mm')
view(73, 12)
colorbar

N = gen.nxGasSupply;

X = reshape(model.(gs).grid.cells.centroids(:, 1), N, [])/(milli*meter);
Y = reshape(model.(gs).grid.cells.centroids(:, 2), N, [])/(milli*meter);

figure('position', [1290, 755, 1275, 559])

val = state.(gs).massfractions{1};
Z = reshape(val, N, []);

surf(X, Y, Z, 'edgecolor', 'none');
colorbar
title('Mass Fraction H2O');
xlabel('x / mm')
ylabel('y / mm')
view([50, 51]);

figure('position', [1290, 755, 1275, 559])

val = state.(gs).massfractions{2};
Z = reshape(val, N, []);

surf(X, Y, Z, 'edgecolor', 'none');
colorbar
title('Mass Fraction O2');
xlabel('x / mm')
ylabel('y / mm')
view([50, 51]);

figure('position', [1290, 755, 1275, 559])

val = state.(gs).pressure;
Z = reshape(val, N, []);

surf(X, Y, Z/barsa, 'edgecolor', 'none');
colorbar
title('Pressure / bar');
xlabel('x / mm')
ylabel('y / mm')
view([50, 51]);

figure('position', [1290, 755, 1275, 559])

val = state.(gs).density;
Z = reshape(val, N, []);

surf(X, Y, Z, 'edgecolor', 'none');
colorbar
title('Density kg/m^3');
xlabel('x / mm')
ylabel('y / mm')
view([50, 51]);


% Current in anode

i = state.Electrolyser.Anode.i;

ind   = model.Electrolyser.couplingTerms{1}.couplingfaces(:, 2);
yc    = model.Electrolyser.Electrolyte.grid.faces.centroids(ind, 2);
areas = model.Electrolyser.Electrolyte.grid.faces.areas(ind);

u = ampere/((centi*meter)^2);
i = (i./areas)/u;

figure
plot(yc/(milli*meter), i);
title('Current in Anode / A/cm^2')
xlabel('height / mm')


% Current in anode

iHp = state.Electrolyser.Anode.iHp;

ind   = model.Electrolyser.couplingTerms{1}.couplingfaces(:, 2);
yc    = model.Electrolyser.Electrolyte.grid.faces.centroids(ind, 2);
areas = model.Electrolyser.Electrolyte.grid.faces.areas(ind);

u = ampere/((centi*meter)^2);
iHp = (iHp./areas)/u;

figure
plot(yc/(milli*meter), iHp);
title('iHp in Anode / A/cm^2')
xlabel('height / mm')

% Faradic effect in Anode

drivingForces.src = @(time) controlfunc(time, I, timeswitch, totaltime, 'order', 'I-first');
state = model.evalVarName(state, 'Electrolyser.Anode.iHp', {{'drivingForces', drivingForces}});

i   = state.Electrolyser.Anode.i;
iHp = state.Electrolyser.Anode.iHp;

ind = model.Electrolyser.couplingTerms{1}.couplingfaces(:, 2);
yc  = model.Electrolyser.Electrolyte.grid.faces.centroids(ind, 2);

figure
plot(yc/(milli*meter), iHp./i);
title('Faradic efficiency')
xlabel('height / mm')


%%

figure
plotToolbar(model.GasSupply.grid, states);

%%
