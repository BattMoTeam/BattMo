function [voltage, U1, U2] = simulateGittEcm(parameters, time, current, soc, ocv)
%SIMULATEGITTECM Simulate a 2-RC ECM with exact discrete state updates.
%
% Current follows BattMo's convention: positive current discharges the cell.

    time = double(time(:));
    current = double(current(:));
    soc = double(soc(:));
    ocv = double(ocv(:));
    n = numel(time);
    assert(numel(current) == n && numel(soc) == n && numel(ocv) == n, ...
           'Time, current, SOC and OCV must have equal lengths.');
    assert(all(diff(time) >= 0), 'Time must be nondecreasing.');

    % Seed the state arrays from a parameter so the same implementation
    % works for both doubles and MRST ADI variables.
    seed = 0*parameters.R0;
    U1 = seed + zeros(n, 1);
    U2 = seed + zeros(n, 1);
    U1(1) = parameters.U10;
    U2(1) = parameters.U20;
    socReference = mean(soc);
    voltage = seed + zeros(n, 1);

    for k = 1:n
        voltage(k) = ocv(k) + parameters.ocvOffset + ...
                     parameters.ocvSlope*(soc(k) - socReference) - ...
                     parameters.R0*current(k) - U1(k) - U2(k);
        if k < n
            dt = time(k + 1) - time(k);
            alpha1 = exp(-dt./parameters.tau1);
            alpha2 = exp(-dt./parameters.tau2);
            U1(k + 1) = alpha1*U1(k) + ...
                        parameters.R1*(1 - alpha1)*current(k);
            U2(k + 1) = alpha2*U2(k) + ...
                        parameters.R2*(1 - alpha2)*current(k);
        end
    end
end
