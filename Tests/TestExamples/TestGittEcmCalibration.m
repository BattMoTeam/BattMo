classdef TestGittEcmCalibration < matlab.unittest.TestCase

    properties
        dataFile
    end

    methods (TestClassSetup)

        function setupBattMo(test)
            repositoryDirectory = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            run(fullfile(repositoryDirectory, 'startupBattMo.m'));
            test.dataFile = fullfile(repositoryDirectory, 'Examples', 'Advanced', ...
                'Impedance', 'data', 'gitt', 'gitt_measurements.csv');
        end

    end

    methods (Test)

        function testNonUniformTimeWeights(test)
            time = [0; 1; 3; 10];
            weights = gittTimeWeights(time);
            test.verifyEqual(weights, [0.05; 0.15; 0.45; 0.35], ...
                             'AbsTol', 1e-14);
            test.verifyEqual(sum(weights), 1, 'AbsTol', 1e-14);
            residual = [0; 0; 1; 1];
            test.verifyEqual(sum(weights.*residual.^2), 0.8, ...
                             'AbsTol', 1e-14);

            groups = ["Pulse"; "Pulse"; "Rest"; "Rest"];
            groupedWeights = gittTimeWeights(time, groups);
            test.verifyEqual(sum(groupedWeights(groups == "Pulse")), 0.5, ...
                             'AbsTol', 1e-14);
            test.verifyEqual(sum(groupedWeights(groups == "Rest")), 0.5, ...
                             'AbsTol', 1e-14);
        end

        function testExactVariableStepSimulationAndAd(test)
            raw = [0.02; 0.01; 5; 0.03; 50; 0; 0; 0; 0];
            parametersAd = initVariablesADI(raw);
            parameters = parameterStruct(parametersAd);
            time = [0; 0.1; 1.3; 4.8; 10];
            current = [0; 1; 1; 0; 0];
            soc = 0.5*ones(size(time));
            ocv = 3.7*ones(size(time));

            voltageAd = simulateGittEcm(parameters, time, current, soc, ocv);
            expected = simulateManually(parameterStruct(raw), time, current, ocv);

            test.verifyClass(voltageAd, 'ADI');
            test.verifyEqual(voltageAd.val, expected, 'AbsTol', 1e-13);
            test.verifySize(voltageAd.jac{1}, [numel(time), numel(raw)]);
            test.verifyEqual(full(voltageAd.jac{1}(:, 1)), -current, ...
                             'AbsTol', 1e-13);
        end

        function testMeasuredCsvCalibration(test)
            data = loadGittData(test.dataFile);
            test.verifyEqual(height(data), 560779);
            test.verifyGreaterThan(max(diff(data.Time_s))/min(diff(data.Time_s)), 1e3);

            segments = segmentGittData(data);
            test.verifyEqual(numel(segments), 483);
            test.verifyTrue(any(strcmp({segments.direction}, 'Charge')));
            test.verifyTrue(any(strcmp({segments.direction}, 'Discharge')));
            charge = segments(find(strcmp({segments.direction}, 'Charge'), 1));
            discharge = segments(find(strcmp({segments.direction}, 'Discharge'), 1));
            test.verifyLessThan(median(charge.current(charge.pulseIndices)), 0);
            test.verifyGreaterThan(median(discharge.current(discharge.pulseIndices)), 0);

            options = struct('socTargets', [0.3, 0.7], ...
                             'maximumOcvIterations', 1, ...
                             'maximumOptimizerIterations', 3, ...
                             'multiStarts', 1, ...
                             'maximumFitSamples', 80);
            result = calibrateEcmFromGitt(data, options);

            test.verifyEqual(height(result.finalFits), 4);
            test.verifyTrue(all(result.finalFits.Accepted));
            test.verifyEqual(height(result.validation), 4);
            test.verifyTrue(all(isfinite(result.validation.RMSE2RC_V)));
            test.verifyTrue(all(result.charge.parameterMap.tau2 >= ...
                                2*result.charge.parameterMap.tau1));
            test.verifyTrue(all(result.discharge.parameterMap.tau2 >= ...
                                2*result.discharge.parameterMap.tau1));
            test.verifyTrue(isfield(result.charge.jsonstruct, 'R2'));
            test.verifyTrue(isfield(result.discharge.jsonstruct, 'C2'));
            test.verifyEqual(result.optimizationMethod, ...
                             '2-RC: unitBoxBFGS with MRST ADI gradients');

            inputParameters = EquivalentCircuitModelInputParams( ...
                result.charge.jsonstruct);
            model = EquivalentCircuitModel(inputParameters);
            test.verifyGreaterThan(model.R2func(0.5), 0);
            test.verifyGreaterThan(model.C2func(0.5), 0);
        end

    end

end

function parameters = parameterStruct(values)
    parameters = struct('R0', values(1), 'R1', values(2), ...
                        'tau1', values(3), 'R2', values(4), ...
                        'tau2', values(5), 'U10', values(6), ...
                        'U20', values(7), 'ocvOffset', values(8), ...
                        'ocvSlope', values(9));
end

function voltage = simulateManually(parameters, time, current, ocv)
    U1 = parameters.U10;
    U2 = parameters.U20;
    voltage = zeros(size(time));
    for k = 1:numel(time)
        voltage(k) = ocv(k) - parameters.R0*current(k) - U1 - U2;
        if k < numel(time)
            dt = time(k + 1) - time(k);
            alpha1 = exp(-dt/parameters.tau1);
            alpha2 = exp(-dt/parameters.tau2);
            U1 = alpha1*U1 + parameters.R1*(1 - alpha1)*current(k);
            U2 = alpha2*U2 + parameters.R2*(1 - alpha2)*current(k);
        end
    end
end
