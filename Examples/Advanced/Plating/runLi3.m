clear
close all

j0_cases = {'Arrhenius', 'Lin'};
styles   = {'-', '--'};
labels   = {'Arrhenius', 'Simplified'};

results = cell(2, 1);

for idx = 1 : numel(j0_cases)

    j0_case = j0_cases{idx};
    
    jsonstruct = parseBattmoJson(fullfile('ParameterData','BatteryCellParameters','LithiumIonBatteryCell','lithium_ion_battery_nmc_graphite.json'));
    jsonstruct_lithium_plating = parseBattmoJson(fullfile('Examples','Advanced','Plating','lithium_plating.json'));

    jsonstruct.use_thermal = false;
    jsonstruct.include_current_collectors = false;

    ne = 'NegativeElectrode'; pe = 'PositiveElectrode';
    co = 'Coating'; am = 'ActiveMaterial'; lp = 'LithiumPlating';

    jsonstruct.(ne).(co).(am).diffusionModelType = 'full';
    jsonstruct.(pe).(co).(am).diffusionModelType = 'full';
    jsonstruct.(ne).(co).(am).useLithiumPlating = false;
    jsonstruct.(ne).(co).(am).isRootSimulationModel = true;

    jsonstruct_lithium_plating.LithiumPlating.kInter = 6e3;
    jsonstruct.(ne).(co).(am).LithiumPlating = jsonstruct_lithium_plating.LithiumPlating;

    jsonstruct.Control.controlPolicy = 'CCCharge';
    
    switch j0_case
      case 'Arrhenius'
        % default case, do nothing
      case 'Lin'
        
        ecd_func = struct('functionFormat', 'named function', ...
                  'functionName'  , 'exchangeCurrentDensityLin');
        ecd_func.argumentList = {'celyte', 'celde', 'T'};
        
        jsonstruct.(ne).(co).(am).(itf).exchangeCurrentDensity = ecd_func;
      otherwise
        
        error('j0_case not recognized');
        
    end
    
    inputparams = BatteryInputParams(jsonstruct);
    inputparams = inputparams.(ne).(co).(am);
    inputparams.Interface.openCircuitPotential.functionname = 'computeOCP_Graphite_Latz';

    model = ActiveMaterial(inputparams);
    model = model.setupForSimulation();

    model.Interface.j0_case = j0_case;

    sd = 'SolidDiffusion'; itf = 'Interface';
    cElectrolyte = 0.5 * mol / litre;
    phiElectrolyte = 0;
    T = 298;

    cInit = model.(itf).guestStoichiometry0 * model.(itf).saturationConcentration;
    N = model.(sd).N;
    initState.(sd).c = cInit * ones(N, 1);
    initState.(sd).cSurface = cInit;

    initState.T = T;
    initState.(itf).cElectrolyte = cElectrolyte;
    initState.(itf).phiElectrolyte = phiElectrolyte;
    initState = model.evalVarName(initState, {itf, 'OCP'});
    OCP = initState.(itf).OCP;
    initState.E = OCP + phiElectrolyte;

    F = model.(itf).constants.F;
    R = model.(itf).constants.R;

    if model.useLithiumPlating
        nPl0 = model.(lp).nPl0;
        r = model.(lp).particleRadius;
        vf = model.(lp).volumeFraction;
        platedConcentration0 = nPl0 * vf / ((4/3)*pi*r^3);
        platedConcentrationInit = platedConcentration0 / (exp((F*OCP)/(R*T)) - 1)^(1/4);
        model.(lp).platedConcentrationRef = platedConcentrationInit;

        initState.(lp).platedConcentrationNorm = platedConcentrationInit / model.(lp).platedConcentrationRef;
        initState.(lp).phiSolid = initState.E;
        initState.(lp).phiElectrolyte = phiElectrolyte;
        initState.(lp).cElectrolyte = cElectrolyte;
        initState.(lp).nSEI = 0;
    end

    Iref = 5e-12;
    Imax = 50 * Iref;
    total = 1*hour*(Iref/Imax);
    dt = total / 100;

    control.src = @(t) rampupControl(t, 1*second*(Iref/Imax), -Imax);
    cmax = model.(itf).guestStoichiometry100 * model.(itf).saturationConcentration;
    control.stopFunction = @(model, state, ~) (state.(sd).cSurface >= cmax);
    schedule = struct('control', control, 'step', struct('val', dt*ones(100, 1), 'control', ones(100, 1)));

    scalingparams = struct('I', Imax, 'elyteConcentration', cElectrolyte);
    if model.useLithiumPlating
        scalingparams.platedConcentration = platedConcentrationInit;
    end
    model = model.setupScalings(scalingparams);

    nls = NonLinearSolver();
    nls.errorOnFailure = false;
    model.nonlinearTolerance = 1e-2;
    model.verbose = false;

    [~, states, ~] = simulateScheduleAD(initState, model, schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls);

    states = states(~cellfun(@isempty, states));
    time = cellfun(@(s) s.time, states);

    varsToEval = { ...
        {'Interface', 'eta'}, ...
        {'Interface', 'intercalationFlux'}, ...
        {'Interface', 'j0'}, ...
        {'SolidDiffusion', 'cSurface'}};

    for k = 1:numel(states)
        for var = 1:numel(varsToEval)
            states{k} = model.evalVarName(states{k}, varsToEval{var});
        end
    end

    for v = 1:length(varsToEval)
        field = varsToEval{v}{end};
        results{idx}.(field) = cellfun(@(s) getfield(s.(varsToEval{v}{1}), field), states);
    end
    results{idx}.time = time;
end

varnames = {'eta', 'intercalationFlux', 'j0', 'cSurface'};

ylabels = {'\eta [V]', 'interflux', 'j0', 'c_{surface} [mol/L]'};

for i = 1:length(varnames)
    figure;
    hold on;
    for idx = 1:2
        plot(results{idx}.time/hour, results{idx}.(varnames{i}), ...
             'LineStyle', styles{idx}, 'LineWidth', 2, ...
             'DisplayName', labels{idx});
    end
    xlabel('Time [hour]');
    ylabel(ylabels{i});
    title(varnames{i});
    legend show;
    grid on;
    hold off;
end
