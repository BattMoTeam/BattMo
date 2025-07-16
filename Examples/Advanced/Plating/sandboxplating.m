%% run stand-alone active material model with lithium plating

clear all
close all

%% Setup the properties of Li-ion battery materials and cell design

jsonstruct = parseBattmoJson(fullfile('ParameterData','BatteryCellParameters','LithiumIonBatteryCell','lithium_ion_battery_nmc_graphite.json'));

inputparams = BatteryInputParams(jsonstruct);

use_cccv = true;
if use_cccv
    inputparams.SOC = 0;
    cccvstruct = struct( 'controlPolicy'     , 'CCCV'       , ...
                         'initialControl'    , 'charging', ...
                         'numberOfCycles'    , 2            , ...
                         'CRate'             , 1.5          , ...
                         'DRate'             , 1            , ...
                         'lowerCutoffVoltage', 3            , ...
                         'upperCutoffVoltage', 4            , ...
                         'dIdtLimit'         , 1e-2         , ...
                         'dEdtLimit'         , 1e-4);
    cccvinputparams = CcCvControlModelInputParams(cccvstruct);
    inputparams.Control = cccvinputparams;
end

ne      = 'NegativeElectrode';
pe      = 'PositiveElectrode';
lp      = 'LithiumPlating';
elyte   = 'Electrolyte';
thermal = 'ThermalModel';
co      = 'Coating';
am      = 'ActiveMaterial';
itf     = 'Interface';
sd      = 'SolidDiffusion';
ctrl    = 'Control';
cc      = 'CurrentCollector';

jsonstruct.use_thermal = false;
jsonstruct.include_current_collectors = false;


jsonstruct.(ne).(co).(am).diffusionModelType = 'full';
jsonstruct.(pe).(co).(am).diffusionModelType = 'full';

jsonstruct.(ne).(co).(am).useLithiumPlating = false;

% Flag pour modèle stand-alone
jsonstruct.(ne).(co).(am).isRootSimulationModel = true;

jsonstruct_lithium_plating = parseBattmoJson(fullfile('Examples', 'Advanced', 'Plating', 'lithium_plating.json'));

jsonstruct.(ne).(co).(am).LithiumPlating = jsonstruct_lithium_plating.LithiumPlating;

% jsonstruct.(ne).(co).(am).(sd).referenceDiffusionCoefficient = 1e-10;

scenario = 'charge';

%% following is not used at particle level (but necessary to initiliaze full battery below)
switch scenario
  case 'charge'
    jsonstruct.Control.controlPolicy = 'CCCharge';
  case 'discharge'
    jsonstruct.Control.controlPolicy = 'CCDischarge';
  otherwise
    error('scenario not recognized');
end

% Setup InputParams
inputparams = BatteryInputParams(jsonstruct);
inputparams = inputparams.(ne).(co).(am);
inputparams.Interface.openCircuitPotential.functionname = 'computeOCP_Graphite_Latz';

%% Setup the model

model = ActiveMaterial(inputparams);

%% Equip model for simulation

model = model.setupForSimulation();

%% Setup initial state

sd  = 'SolidDiffusion';
itf = 'Interface';

cElectrolyte   = 5e-1*mol/litre;
phiElectrolyte = 0;
T              = 50;

switch scenario
  case 'charge'
    % cElectrodeInit = (model.(itf).guestStoichiometry0)*(model.(itf).saturationConcentration);
    cElectrodeInit = 25*mol/litre;
  case 'discharge'
    cElectrodeInit = (model.(itf).guestStoichiometry100)*(model.(itf).saturationConcentration);
  otherwise
    error('scenario not recognized');
end

N = model.(sd).N;
initState.(sd).c        = cElectrodeInit*ones(N, 1);
initState.(sd).cSurface = cElectrodeInit;

initState.T = T;
initState.(itf).cElectrolyte   = cElectrolyte;
initState.(itf).phiElectrolyte = phiElectrolyte;

initState = model.evalVarName(initState, {itf, 'OCP'});
OCP = initState.(itf).OCP;
initState.E = OCP + phiElectrolyte;


F = model.(itf).constants.F;
R = model.(itf).constants.R;

if model.useLithiumPlating
    nPl0                 = model.LithiumPlating.nPl0;
    r                    = model.LithiumPlating.particleRadius;
    vf                   = model.LithiumPlating.volumeFraction;
    platedConcentration0 = nPl0 * vf / ((4/3)*pi*r^3);
    
    %initialisation so that overpotential are =0 at the beginning
    platedConcentrationInit = platedConcentration0/(exp((F*OCP)/(R*T)) - 1)^(1/4);

    model.(lp).platedConcentrationRef = platedConcentrationInit;

    initState.(lp).platedConcentrationNorm = platedConcentrationInit / model.(lp).platedConcentrationRef;
    initState.(lp).phiSolid            = initState.E;
    initState.(lp).phiElectrolyte      = phiElectrolyte;
    initState.(lp).cElectrolyte        = cElectrolyte;
    initState.(lp).nSEI                = 0;
end


%% setup non-linear solver

nls = NonLinearSolver();
nls.errorOnFailure  = false;
nls.maxTimestepCuts = 20;

model.nonlinearTolerance = 1e-6;

% --- Paramètres de simulation ---
Iref_vals = 1e-5*[1e-12, 3e-12, 5e-12];  % différentes intensités à tester
varnames = {'eta', 'etaPlating', 'etaChemical', ...
            'platingFlux', 'chemicalFlux', 'intercalationFlux', ...
            'surfaceCoverage', 'platedConcentration'};

allResults = struct();

% --- Boucle sur les intensités ---
% Iref_vals = Iref_vals(1);

for i = 1:length(Iref_vals)
    
    Iref = Iref_vals(i);
    Imax = Iref;

    total = 3e-5*hour*(Iref/Imax)*1e9;
    n     = 500;
    dt    = total/n;
    step  = struct('val', dt*ones(n, 1), 'control', ones(n, 1));

    tup = 1e-1*total;

    switch scenario
        case 'charge'
            srcfunc = @(time) rampupControl(time, tup, -Imax); %0 pour tourner à vide
            % srcfunc = @(time) 0; %0 pour tourner à vide
            % srcfunc = @(time) -Imax;
            cmax = (model.(itf).guestStoichiometry100)*(model.(itf).saturationConcentration);
            % control.stopFunction = @(model, state, state0_inner) (state.(sd).cSurface >= cmax);
        case 'discharge'
            srcfunc = @(time) rampupControl(time, tup, Imax); %0 pour tourner à vide
            cmin = (model.(itf).guestStoichiometry0)*(model.(itf).saturationConcentration);
            control.stopFunction = @(model, state, state0_inner) (state.(sd).cSurface <= cmin);
        otherwise
            error('scenario not recognized');
    end

    control.src = srcfunc;

    schedule = struct('control', control, 'step', step);

    scalingparams = struct('I'                  , Imax                             , ...
                           'elyteConcentration' , initState.(itf).cElectrolyte);

    if model.useLithiumPlating
        scalingparams.platedConcentration = platedConcentrationInit;
    end


    model = model.setupScalings(scalingparams);
    % Rescale
    scalingparams.I = Imax;
    if model.useLithiumPlating
        scalingparams.platedConcentration = platedConcentrationInit;
    end
    model = model.setupScalings(scalingparams);

    % Re-initialiser état
    initStateCurrent = initState;

    % Simulation
    [~, states, ~] = simulateScheduleAD(initStateCurrent, model, schedule, ...
                                        'OutputMinisteps', true, 'NonLinearSolver', nls);

    % Nettoyage
    states = states(cellfun(@(s) ~isempty(s), states));
    states = cellfun(@(s) model.evalVarName(s, {sd, 'cAverage'}), states, 'uniformoutput', false);
    
    time = cellfun(@(s) s.time, states);
    time = [0; time];
    
    Q = -diff(time).* srcfunc(time(2 : end));
    Q = cumsum(Q);
    
    Q_mAh = Q / 3.6;            % milliampère-heure

    % % Évaluer les variables d’intérêt
    % varsToEval = {{'Interface'     , 'eta'}         , ...
    %               {'LithiumPlating', 'etaPlating'}  , ...
    %               {'LithiumPlating', 'etaChemical'} , ...
    %               {'Interface'     , 'intercalationFlux'}           , ...
    %               {'LithiumPlating', 'platingFlux'} , ...
    %               {'LithiumPlating', 'chemicalFlux'}, ...
    %               {'LithiumPlating', 'surfaceCoverage'}, ...
    %               {'LithiumPlating', 'platedConcentration'}};
    % 
    % for k = 1:numel(states)
    %     for var = 1:numel(varsToEval)
    %         states{k} = model.evalVarName(states{k}, varsToEval{var});
    %     end
    % end

    % % Extraction des valeurs
    % values = cell(1, numel(varnames));
    % values{1} = cellfun(@(s) s.(itf).eta, states);
    % values{2} = cellfun(@(s) s.(lp).etaPlating, states);
    % values{3} = cellfun(@(s) s.(lp).etaChemical, states);
    % values{4} = cellfun(@(s) s.(lp).platingFlux, states);
    % values{5} = cellfun(@(s) s.(lp).chemicalFlux, states);
    % values{6} = cellfun(@(s) s.(itf).intercalationFlux, states);
    % values{7} = cellfun(@(s) s.(lp).surfaceCoverage, states);
    % values{8} = cellfun(@(s) s.(lp).platedConcentration, states);

    cSurface = cellfun(@(s) s.(sd).cSurface, states);
    cAverage = cellfun(@(s) s.(sd).cAverage, states);
    allResults(i).cSurface = cSurface;
    allResults(i).cAverage = cAverage;
    % Stockage
    allResults(i).Iref = Iref;
    allResults(i).Q_mAh = Q_mAh;
    % for j = 1:numel(varnames)
    %     allResults(i).(varnames{j}) = values{j};
    % end
end

linestyles = {'-', '--', ':', '-.'};        
colors     = lines(length(Iref_vals));       

% for j = 1:numel(varnames)
%     figure
%     hold on
%     for i = 1:length(Iref_vals)
%         ls = linestyles{mod(i-1, length(linestyles)) + 1};
%         plot(allResults(i).Q_mAh, allResults(i).(varnames{j}), ...
%             'LineStyle', ls, ...
%             'Color', colors(i,:), ...
%             'DisplayName', sprintf('I_{ref} = %.0e A', allResults(i).Iref));
%     end
%     xlabel('Charge passed [mAh]');
%     ylabel(varnames{j});
%     title([varnames{j} ' vs Charge passed']);
%     legend show
%     grid on
% end

figure
hold on
for i = 1:length(Iref_vals)
    ls = linestyles{mod(i-1, length(linestyles)) + 1};
    plot(allResults(i).Q_mAh, allResults(i).cSurface / (mol/litre), ...
        'LineStyle', ls, ...
        'Color', colors(i, :), ...
        'DisplayName', sprintf('I_{ref} = %.0e A', allResults(i).Iref));
end
xlabel('Charge passed [mAh]');
ylabel('Surface concentration [mol/L]');
title('Surface concentration vs Charge passed');
legend show
grid on
