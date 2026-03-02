%% Multi-reactionRatePlating and reactionRateChemicalIntercalation simulation with lithium plating

clear
close all

%% Define reactionRatePlating and reactionRateChemicalIntercalation values to test
reactionRatePlating_values = logspace(-11, -9, 10);
reactionRateChemicalIntercalation_values = logspace(-15, -9, 40);
n_reactionRatePlating = length(reactionRatePlating_values);
n_reactionRateChemicalIntercalation = length(reactionRateChemicalIntercalation_values);

results = cell(n_reactionRatePlating, n_reactionRateChemicalIntercalation);
line_styles = {'-', '--', ':'};

for i = 1:n_reactionRatePlating
    for j = 1:n_reactionRateChemicalIntercalation
        reactionRatePlating = reactionRatePlating_values(i);
        reactionRateChemicalIntercalation = reactionRateChemicalIntercalation_values(j);
        fprintf('\nRunning simulation (%d/%d, %d/%d): reactionRatePlating = %.1e, reactionRateChemicalIntercalation = %.1e\n', ...
            i, n_reactionRatePlating, j, n_reactionRateChemicalIntercalation, reactionRatePlating, reactionRateChemicalIntercalation);

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

        jsonstruct_lithium_plating.LithiumPlating.reactionRatePlating = reactionRatePlating;
        jsonstruct_lithium_plating.LithiumPlating.reactionRateChemicalIntercalation = reactionRateChemicalIntercalation;
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
            thresholdParameter = model.(lp).thresholdParameter;
            r = model.(lp).particleRadius;
            vf = model.(lp).volumeFraction;
            platedConcentration0 = thresholdParameter * vf / ((4/3)*pi*r^3);
            platedConcentrationInit = platedConcentration0 / (exp((F*OCP)/(R*T)) - 1)^(1/4);
            model.(lp).platedReferenceConcentration = platedConcentrationInit;

            initState.(lp).platedConcentrationNorm = platedConcentrationInit / model.(lp).platedReferenceConcentration;
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
        results{i,j}.reactionRatePlating = reactionRatePlating;
        results{i,j}.reactionRateChemicalIntercalation = reactionRateChemicalIntercalation;
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
    for i = 1:n_reactionRatePlating
        for j = 1:n_reactionRateChemicalIntercalation
            plot(results{i,j}.time/hour, results{i,j}.(labels{i_label}), ...
                'LineStyle', '-', ...
                'LineWidth', 1.5, ...
                'DisplayName', sprintf('reactionRatePlating=%.1e, reactionRateChemicalIntercalation=%.1e', results{i,j}.reactionRatePlating, results{i,j}.reactionRateChemicalIntercalation));
        end
    end
    xlabel('Time [hour]');
    ylabel(ylabels{i_label});
    title(labels{i_label});
    legend show;
    grid on;
    hold off;
end
