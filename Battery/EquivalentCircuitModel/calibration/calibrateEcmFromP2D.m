function jsonstruct_ecm = calibrateEcmFromP2D(jsonstruct_p2d, soc_values)
    % CALIBRATEECMFROMP2D Calibrates a 2-RC Equivalent Circuit Model (ECM)
    % from a physics-based P2D model across a range of State of Charge (SOC) points.
    %
    % Inputs:
    %   jsonstruct_p2d - P2D model configuration structure
    %   soc_values           - SOC values
    %
    % Output:
    %   jsonstruct_ecm - Standardized BattMo ECM configuration structure
    
    
    [model, inputparams] = setupModelFromJson(jsonstruct_p2d);
    state_current = model.setupInitialState();
    
    N_steps = length(soc_values);
    
    results_matrix = zeros(N_steps, 6); % Stores [SOC, R0, R1, C1, R2, C2]
    ocv_list = zeros(N_steps, 1);       % Stores Open Circuit Voltage for each SOC
    
    freq = logspace(-4, 2, 50)./(2*pi);
    
    %% 2. Main Loop
    for i = 1 : N_steps
        
        current_soc = soc_values(i);
        fprintf('\n==================================================\n');
        fprintf('  EIS fitting for SOC = %.2f (%d/%d)\n', current_soc, i, N_steps);
        fprintf('==================================================\n');
        
        % Run the Impedance Solver
        options = struct();
        options.stateInitialization.initializationSetup = 'soc';
        options.stateInitialization.soc                 = current_soc;
        
        impsolv = ImpedanceSolver(inputparams, options);
        ocv_list(i) = impsolv.getOCP();
        
        fprintf('Linearization...\n');
        Z = impsolv.computeImpedance(freq);
        Z_re = real(Z)';
        Z_im = imag(Z)';
        
        omega = (2*pi).*freq;
        
        % Run parameter identification
        eisdata = struct('omega', omega, ...
                         'Z_re_exp', Z_re, ...
                         'Z_im_exp', Z_im);

        params0 = [0.05052; 1.12673; 59119.9; 0.03155; 11054.0];  % initial condition: C1>2*C2
        
        feis = FittingEIS(eisdata, params0);
        [p, ~] = feis.run();

        % Store values in the matrix
        results_matrix(i, 1) = current_soc;
        results_matrix(i, 2:6) = p;
        
    end
    
    %% 3. Post-Processing & Sorting (Ascending SOC Order)
    % BattMo ECM JSON requires OCP and parameters defined with ascending SOC (0 -> 1)
    [~, sort_idx] = sort(soc_values, 'ascend');
    soc_asc = soc_values(sort_idx)';
    ocv_asc = ocv_list(sort_idx)';
    
    results_matrix_sorted = results_matrix(sort_idx, :);
    R0_asc = results_matrix_sorted(:, 2)';
    R1_asc = results_matrix_sorted(:, 3)';
    C1_asc = results_matrix_sorted(:, 4)';
    R2_asc = results_matrix_sorted(:, 5)';
    C2_asc = results_matrix_sorted(:, 6)';
    
    %% 4. Build jsonstruct_ecm
    jsonstruct_ecm = struct();
    
    % Tabulated Open Circuit Potential (OCP)
    jsonstruct_ecm.OCP.functionFormat = 'tabulated';
    jsonstruct_ecm.OCP.argumentList = {'SOC'};
    jsonstruct_ecm.OCP.dataX = soc_asc;
    jsonstruct_ecm.OCP.dataY = ocv_asc;
    
    % Tabulated R0
    jsonstruct_ecm.R0.functionFormat = 'tabulated';
    jsonstruct_ecm.R0.argumentList = {'SOC'};
    jsonstruct_ecm.R0.dataX = soc_asc;
    jsonstruct_ecm.R0.dataY = R0_asc;
    
    % Tabulated R1 & C1 (First RC branch)
    jsonstruct_ecm.R1.functionFormat = 'tabulated';
    jsonstruct_ecm.R1.argumentList = {'SOC'};
    jsonstruct_ecm.R1.dataX = soc_asc;
    jsonstruct_ecm.R1.dataY = R1_asc;
    
    jsonstruct_ecm.C1.functionFormat = 'tabulated';
    jsonstruct_ecm.C1.argumentList = {'SOC'};
    jsonstruct_ecm.C1.dataX = soc_asc;
    jsonstruct_ecm.C1.dataY = C1_asc;
    
    % Tabulated R2 & C2 (Second RC branch)
    jsonstruct_ecm.R2.functionFormat = 'tabulated';
    jsonstruct_ecm.R2.argumentList = {'SOC'};
    jsonstruct_ecm.R2.dataX = soc_asc;
    jsonstruct_ecm.R2.dataY = R2_asc;
    
    jsonstruct_ecm.C2.functionFormat = 'tabulated';
    jsonstruct_ecm.C2.argumentList = {'SOC'};
    jsonstruct_ecm.C2.dataX = soc_asc;
    jsonstruct_ecm.C2.dataY = C2_asc;
    
    % Print and save summary
    ecm_table = array2table(results_matrix_sorted, 'VariableNames', ...
        {'SOC', 'R0', 'R1', 'C1', 'R2', 'C2'});
    disp(ecm_table);
    
end

%{
Copyright 2021-2026 SINTEF Industry, Sustainable Energy Technology
and SINTEF Digital, Mathematics & Cybernetics.

This file is part of The Battery Modeling Toolbox BattMo

BattMo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BattMo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BattMo.  If not, see <http://www.gnu.org/licenses/>.
%}
