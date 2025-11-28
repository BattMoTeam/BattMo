%% Multi-kPl and kChInt simulation with lithium plating

clear
close all

%% Define kPl and kChInt values to test
kPl_values = logspace(-11, -9, 10);
kChInt_values = logspace(-15, -9, 40);
n_kPl = length(kPl_values);
n_kChInt = length(kChInt_values);

results = cell(n_kPl, n_kChInt);
line_styles = {'-', '--', ':'};

for i = 1:n_kPl
    for j = 1:n_kChInt
        kPl = kPl_values(i);
        kChInt = kChInt_values(j);
        fprintf('\nRunning simulation (%d/%d, %d/%d): kPl = %.1e, kChInt = %.1e\n', ...
            i, n_kPl, j, n_kChInt, kPl, kChInt);

        jsonstruct = parseBattmoJson(fullfile('ParameterData','BatteryCellParameters','LithiumIonBatteryCell','lithium_ion_battery_nmc_graphite.json'));
        jsonstruct_lithium_plating = parseBattmoJson(fullfile('Examples','Advanced','Plating','lithium_plating.json'));

        jsonstruct.use_thermal = false;
        jsonstruct.include_current_collectors = false;

        ne = 'NegativeElectrode'; pe = 'PositiveElectrode';
        co = 'Coating'; am = 'ActiveMaterial'; lp = 'LithiumPlating';

        jsonstruct.(ne).(co).(am).diffusionModelType = 'full';
        jsonstruct.(pe).(co).(am).diffusionModelType = 'full';
        jsonstruct.(ne).(co).(am).useLithiumPlating = true;
        jsonstruct.(ne).(co).(am).isRootSimulationModel = true;

        jsonstruct_lithium_plating.LithiumPlating.kPl = kPl;
        jsonstruct_lithium_plating.LithiumPlating.kChInt = kChInt;
        jsonstruct.(ne).(co).(am).LithiumPlating = jsonstruct_lithium_plating.LithiumPlating;

        jsonstruct.Control.controlPolicy = 'CCCharge';

        inputparams = BatteryInputParams(jsonstruct);
        inputparams = inputparams.(ne).(co).(am);
        inputparams.Interface.openCircuitPotential.functionname = 'computeOCP_Graphite_Latz';

        model = ActiveMaterial(inputparams);
        model = model.setupForSimulation();

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

        varsToEval = {
            {'SolidDiffusion', 'cSurface'}, ...
            {'Interface', 'eta'}, ...
            {'Interface', 'intercalationFlux'}, ...
            {'LithiumPlating', 'etaPlating'}, ...
            {'LithiumPlating', 'etaChemical'}, ...
            {'LithiumPlating', 'platingFlux'}, ...
            {'LithiumPlating', 'chemicalFlux'}, ...
            {'LithiumPlating', 'surfaceCoverage'}, ...
            {'LithiumPlating', 'platedConcentrationNorm'}};

        for k = 1:numel(states)
            for var = 1:numel(varsToEval)
                states{k} = model.evalVarName(states{k}, varsToEval{var});
            end
        end

        for v = 1:length(varsToEval)
            varname = varsToEval{v}{end};
            results{i,j}.(varname) = cellfun(@(s) getfield(s.(varsToEval{v}{1}), varname), states);
        end
        results{i,j}.time = time;
        results{i,j}.kPl = kPl;
        results{i,j}.kChInt = kChInt;
    end
end

%% Optional: Example plot (you can customize further)
labels = {'cSurface', 'eta', 'etaPlating', 'etaChemical', 'platingFlux', ...
          'chemicalFlux', 'intercalationFlux', 'surfaceCoverage', 'platedConcentrationNorm'};
ylabels = {'Surface concentration [mol/m3]', 'Overpotential [V]', 'Plating η [V]', 'Chemical η [V]', ...
           'Plating Flux [mol/m²/s]', 'Chemical Flux [mol/m²/s]', ...
           'Intercalation Flux [mol/m²/s]', 'Surface Coverage', ...
           'Normalized Plated Conc.'};

for i_label = 1:length(labels)
    figure;
    hold on;
    for i = 1:n_kPl
        for j = 1:n_kChInt
            plot(results{i,j}.time/hour, results{i,j}.(labels{i_label}), ...
                'LineStyle', '-', ...
                'LineWidth', 1.5, ...
                'DisplayName', sprintf('kPl=%.1e, kChInt=%.1e', results{i,j}.kPl, results{i,j}.kChInt));
        end
    end
    xlabel('Time [hour]');
    ylabel(ylabels{i_label});
    title(labels{i_label});
    legend show;
    grid on;
    hold off;
end
